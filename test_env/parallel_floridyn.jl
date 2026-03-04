using Timers
using Evolutionary
using OhMyThreads: TaskLocalValue
using Infiltrator
using Profile
using BenchmarkTools
using LinearAlgebra
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end
println()
println("AAAAA note: launching dev branch")
println("AAAAA note: launching dev branch")
println("AAAAA note: launching dev branch")
println()
#get the visualisation and settings
#_ , vis_file = get_default_project()[2:3]
vis_file = "data/vis_default.yaml"
settings_file = "data/REALWF_CONTROLTEST_VEER.yaml"   #custom data file with veer specification
vis = Vis(vis_file)
plt=nothing

include("../examples/remote_plotting.jl")
include("optimisationstructs.jl")
include("functions.jl")


# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
con.yaw_data = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate yaw matrix with zeros to prevent error on 'nothing'
#run initial conditions
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 


#thread-safety for parallel execution     
if !@isdefined(tlv)
    const tlv = TaskLocalValue{NamedTuple}() do
        (
            plt=deepcopy(plt), set=deepcopy(set), wf=deepcopy(wf),
            wind=deepcopy(wind), sim=deepcopy(sim), floris=deepcopy(floris),
            floridyn=deepcopy(floridyn), vis=deepcopy(vis), con=deepcopy(con),

            power_vector = zeros(wf.nT * sim.n_sim_steps), #preallocate buffer to not return DataFrame structures
            yaws_with_time_buffer = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1), #preallocate buffer for yaw matrix to avoid allocations in construct_yaw_matrix_dynamic
            wind_dirs_buffer = zeros(sim.end_time - sim.start_time + 1) #preallocate buffer for wind directions at each time step to avoid allocations in get_wind_at_t
        )
    end
end


#------------------------------------------------------------------------------------------------------------
#OPTIMISATION PART


#optimisation constraints
set_num_yaw_changes = 4  #N
set_num_optimiser_runs = 2 #number of automatic restarts for CMA-ES
set_max_yaw_rate = 1.0 #deg/s
set_sigma0 =  0.05         # 0.05 for time optim, 0.01 works well in second run for yaws!!     # set to 30% of the search range, and for yaw convergence: first 0.1 for time , then 0.03 for yaws
set_max_yaw_misalignment = 45.0 #deg, for penalising large yaw angles in the cost function, for stability and convergence reasons
set_lambda_l1 = 0.0 #1e3  #units: cost PER DEGREE, per turbine, PER SECOND #relative to the beneficial term average kW per turbine #typical value 1e3, can play around with this
set_lambda_l1_hard_limit = Inf #a limit on the maximum total yaw change in a simulation, in degrees
set_objective = totalEnergyObjective      #totalEnergyObjective or powerTrackingObjective 
#set_individual_turbine_switching = false  #whether to allow individual turbine switching or only simultaneous switching, for simplicity we will stick with simultaneous switching for now, and can report about this in the thesis
#set_desired_power_curve = ones(sim.n_sim_steps) *  60.0  #find a way to elegantly match these numbers! 



#initial state
x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)   #x0 = result.minimizer  #start from previous result also possible
                        


#hyperparams
set_lambda_multiplier = 2 #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
set_lambda0 = 2 * round(   (4 + 3 * log(wf.nT * (set_num_yaw_changes-1)))      / 2.0   )   #half the offspring 
set_lambda = Int(set_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))




#efficient struct passing
opt_set = OptimisationSettings(
    set_num_yaw_changes, set_num_optimiser_runs, set_max_yaw_rate, set_lambda_multiplier, set_lambda0, set_lambda, set_mu,
    set_sigma0, set_sigma0, set_max_yaw_misalignment, set_lambda_l1, set_lambda_l1_hard_limit, set_objective
)



# set cost function
cost_func = parallel_costfunction(plt, set, wf, wind, sim, con, vis, floridyn, floris, opt_set)


#opts
opts = Evolutionary.Options( 
    iterations = 10,
    abstol = 1e-8,
    reltol = 1e-8,
    show_trace = true,
    show_every = 5, 
    store_trace = true,
    parallelization = :thread  #serial   #thread           #multithreading hehe
)


#constraints
#limit normalised time to [0,1], yaws are free after adding try/catch, a values are bound.
lower_bounds = zeros(set_num_yaw_changes-1)
upper_bounds = ones(set_num_yaw_changes-1)
println("Starting CMA-ES run 1")
println("sigma0=$set_sigma0, lambda=$set_lambda, mu=$set_mu")



begin_time = time()
@time result =  Evolutionary.optimize(cost_func,
                                BoxConstraints(lower_bounds, upper_bounds),
                                x0, 
                                CMAES(  lambda = set_lambda,
                                        mu = set_mu,
                                        sigma0 = set_sigma0),
                                opts)




#secondary runs with other options as well
for run in 2:set_num_optimiser_runs

    #hyperparams
    local_set_lambda_multiplier = set_lambda_multiplier #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
    local_set_lambda0 = set_lambda0 
    local_set_lambda = set_lambda
    local_set_mu = set_mu

    #reduce sigma to stay in local region
    local_set_sigma0 = 0.01

    #print
    println()   
    println("Starting CMA-ES run $run")
    println("sigma0=$local_set_sigma0, lambda=$local_set_lambda, mu=$local_set_mu")

    #set new initial guess and reinitialise cov matrix
    local_x0 = result.minimizer

    @time global result =  Evolutionary.optimize(cost_func,        
                                BoxConstraints(lower_bounds, upper_bounds),
                                local_x0, 
                                CMAES(  lambda = local_set_lambda,
                                        mu = local_set_mu,
                                        sigma0 = local_set_sigma0),
                                opts)
end


end_time = time()
total_time = end_time - begin_time
println()
println("Total optimization time: $(round(total_time, digits=2)) seconds")



#----------------------------------------------------------------------------------

#PLOTTING & POSTPROCESSING







# the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)


#set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();

#AAAAAAAAA note: if we keep a between 0.2 and 0.33, we can approximate with constant lambda of 7.6
#meaning we can do axial induction control in this range without changing the code dynamics!!!
# == TO DO FOR LATER

#run initial conditions
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
wf = initSimulation(wf, sim);
con.yaw_data = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate yaw matrix with zeros to prevent error on 'nothing'


#disable online visualisation
vis.online = false 


#plot flowfield
#construct_yaw_matrix_dynamic!(con.yaw_data, result.minimizer, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
construct_yaw_matrix_dynamic!(con.yaw_data, result.minimizer, sim, wf, opt_set)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)


#include("get_energy.jl")  #for energy tests
#include("controller_output.jl")   #to see and plot controller angles

println()
include("calculate_increase_over_baseline.jl")  #to calculate the increase over baseline for the optimised case, compared to a baseline case with no yawing






#after all is done, perform GC already
GC.gc()



#notes--------------------------------------------------------------------------------

#=






#=   for power maximisation: objective > 0 means failure as problem is not feasible
println()
if result.minimum > 0
    println("ERROR: optimisation did not converge")
end
println("-----------------------")
=#



also some more notes:
20 feb 2026
individual turbine switching has a lower minimum but finding it is too hard for the optimiser,
so we will stick with simultaneous switching for now, and can report about this in the thesis



next steps are minimizing the l1 norm of the yaw actuation on top of the power maximisation. 
the code provides an interface for this via the cost function:
    max_yaw_misalignment is a way of penalising large yaw angles (also for stability)
    


also next time do profiling to test the memory bottleneck with ~17% gc

 
3 march:
added axial induction control! but this type of control needs smoother adjustments of the A matrix as this leads to very spiky power.
to run: 

include parallel_floridyn.jl
include test.jl






=#






