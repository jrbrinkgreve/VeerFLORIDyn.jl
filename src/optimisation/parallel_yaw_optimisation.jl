using Timers
using Evolutionary
using OhMyThreads: TaskLocalValue
using Infiltrator
using Profile
using BenchmarkTools
using LinearAlgebra
using FLORIDyn, TerminalPager, DistributedNext
using Plots

if Threads.nthreads() == 1; using ControlPlots; end


#region fold initialisation and setup
#close old plots:

#get the visualisation and settings
vis_file = "data/vis_default.yaml"
settings_file = "data/REALWF_CONTROLTEST_VEER.yaml"   #   "data/JHTDB_comparison_turbines.yaml"   #"data/REALWF_CONTROLTEST_VEER.yaml"   #custom data file with veer specification
vis = Vis(vis_file)
plt=nothing

include("../../examples/remote_plotting.jl")
include("optimisationstructs.jl")
include("functions.jl")
include("result_generation/helper_functions.jl")


#BLAS/multithreading 
BLAS.set_num_threads(1)

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

if !@isdefined(tlv)
    const tlv = Ref{TaskLocalValue{NamedTuple}}()
end

function init_tlv!()
    tlv[] = TaskLocalValue{NamedTuple}() do
        (
            plt      = deepcopy(plt),
            set      = deepcopy(set),
            wf       = deepcopy(wf),
            wind     = deepcopy(wind),
            sim      = deepcopy(sim),
            floris   = deepcopy(floris),
            floridyn = deepcopy(floridyn),
            vis      = deepcopy(vis),
            con      = deepcopy(con),
            power_vector          = zeros(wf.nT * sim.n_sim_steps),
            yaws_with_time_buffer = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1),
            wind_dirs_buffer      = zeros(sim.end_time - sim.start_time + 1),
        )
    end
end

# (Re)initialise whenever data changes
init_tlv!()





#endregion

#------------------------------------------------------------------------------------------------------------------
#OPTIMISATION PART

#sim 
set_num_yaw_changes = 4        #N
set_max_yaw_misalignment = 25 #deg, for penalising large yaw angles in the cost function, for stability and convergence reasons
set_lambda_l1 = 0            #1e3  #units: cost PER DEGREE, per turbine, PER SECOND #relative to the beneficial term average kW per turbine #typical value 1e3, can play around with this
set_lambda_l1_hard_limit = Inf  #a limit on the maximum total yaw change in a simulation, in degrees
set_max_yaw_rate = 1.0          #deg/s
set_objective = totalEnergyObjective   #totalEnergyObjective or totalEnergyObjectivePlusCurvatureRegulariser
set_num_timesteps_to_skip = 125     #skip the first N timesteps for wake effects to propagate, approx time between wake interactions

#optimiser convergence/ fidelity  
set_num_optimiser_runs = 3          # 4       #number of automatic restarts for CMA-ES
set_cmaes_lambda_multiplier = 3     # 4       #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
set_iterations = 30                 # 50     #number of iterations for CMAES
set_sigma0 = 0.05                   # 0.05   #for time optim, 0.01 works well in second run for yaws!!     # set to 30% of the search range, and for yaw convergence: first 0.05 then 0.01 
set_sigma0_secondary = 0.02         # 0.02   #for second run with yaws, to reduce the search area and converge faster, can play around with this #0.01
set_sigma0_final = 0.01             # 0.01  #for final run 

#output
verbose = true
trace_steps = 5



#initial state
x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)

#other initialisations                                          #change this to rand() .- 0.5
#x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes) + 0.05 * randn(size(generate_initial_guess(sim, wind, wf, set_num_yaw_changes)))   
#x0 = result.minimizer           #start from previous result also possible

#region CMAES prep
#efficient struct passing
opt_set = OptimisationSettings(
    set_num_yaw_changes, set_num_optimiser_runs, set_max_yaw_rate, set_max_yaw_misalignment, set_lambda_l1, set_lambda_l1_hard_limit, set_objective, set_num_timesteps_to_skip
)

#automatic hyperparameter setting

set_lambda0 = 2 * round(   (4 + 3 * log(wf.nT * (set_num_yaw_changes)))      / 2.0   )   #half the offspring 
set_lambda = Int(set_cmaes_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))



#set cost function
cost_func = parallel_costfunction(plt, set, wf, wind, sim, con, vis, floridyn, floris, opt_set)


#opts
opts = Evolutionary.Options( 
    iterations = set_iterations,
    abstol = 1e-2,
    reltol = 1e-5,          #1e-8
    show_trace = verbose,
    show_every = trace_steps, 
    store_trace = true,
    parallelization = :thread  #serial   #thread           #multithreading hehe
)



#constraints
#limit normalised time to [0,1], yaws are free
lower_bounds = zeros(set_num_yaw_changes-1)
upper_bounds = ones(set_num_yaw_changes-1)


GC.gc()


println()
println("Starting CMAES run 1 / $set_num_optimiser_runs")
println("σ_0 = $set_sigma0, λ = $set_lambda, μ = $set_mu")
println()
all_traces = [] #to store traces from multiple runs 
#endregion CMAES prep

#region CMAES calls

begin_time = time()
@time result =  Evolutionary.optimize(cost_func,
                                BoxConstraints(lower_bounds, upper_bounds),
                                x0, 
                                CMAES(  lambda = set_lambda,
                                        mu = set_mu,
                                        sigma0 = set_sigma0),
                                opts)

#store trace
push!(all_traces, result.trace)


#println("Result 1")
#println(result)

#secondary runs with other options as well
for run in 2:set_num_optimiser_runs-1

    #hyperparams
    local local_set_lambda_multiplier = set_cmaes_lambda_multiplier #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
    local local_set_lambda0 = set_lambda0 
    local local_set_lambda = set_lambda
    local local_set_mu = set_mu

    #reduce sigma to stay in local region
    local local_set_sigma0 = set_sigma0_secondary

    #print
    println()   
    println("Starting CMAES run $run / $set_num_optimiser_runs")
    println("σ_0 = $local_set_sigma0, λ = $local_set_lambda, μ = $local_set_mu")
    println()

    #set new initial guess and reinitialise cov matrix
    local_x0 = result.minimizer

    @time global result =  Evolutionary.optimize(cost_func,        
                                BoxConstraints(lower_bounds, upper_bounds),
                                local_x0, 
                                CMAES(  lambda = local_set_lambda,
                                        mu = local_set_mu,
                                        sigma0 = local_set_sigma0),
                                opts)
    push!(all_traces, result.trace) #store trace
    #println("Result $run")
    #println(result)
end



if set_num_optimiser_runs > 1

    #final run with even smaller sigma and possibly different hyperparams as well
    local_set_lambda_multiplier = set_cmaes_lambda_multiplier #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
    local_set_lambda0 = set_lambda0 
    local_set_lambda = set_lambda
    local_set_mu = set_mu

    #reduce sigma to stay in local region
    local_set_sigma0 = set_sigma0_final
    #print
    println()   
    println("Starting CMAES run $set_num_optimiser_runs / $set_num_optimiser_runs")
    println("σ_0 = $set_sigma0_final, λ = $set_lambda, μ = $set_mu")
    println()

    #set new initial guess and reinitialise cov matrix
    final_initialiser = result.minimizer

    @time global result =  Evolutionary.optimize(cost_func,        
                                BoxConstraints(lower_bounds, upper_bounds),
                                final_initialiser, 
                                CMAES(  lambda = set_lambda,
                                        mu = set_mu,
                                        sigma0 = set_sigma0_final),
                                opts)
    push!(all_traces, result.trace) #store trace
end
#println("Result $set_num_optimiser_runs")
#println(result)

end_time = time()
total_time = end_time - begin_time
println()
println("Total optimization time: $(round(total_time, digits=2)) seconds")

#endregion CMAES

#------------------------------------------------------------------------------------------------------------------

#region PLOTTING & POSTPROCESSING



# the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)


#set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();

#run initial conditions
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
wf = initSimulation(wf, sim);
con.yaw_data = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate yaw matrix with zeros to prevent error on 'nothing'


#disable online visualisation
vis.online = false 


val, arg = max_yaw_misalignment(result.minimizer)
println("Maximum yaw misalignment: $(round(val, digits=3)) degrees, at time step $(arg)")


#plot flowfield

#construct_yaw_matrix_dynamic!(con.yaw_data, x0, sim, wf, opt_set);
construct_yaw_matrix_dynamic!(con.yaw_data, result.minimizer, sim, wf, opt_set);
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris);

#top-view
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis);
plot_flow_field(wf, X, Y, Z, vis; msr=AddedTurbulence, plt);



#region cross-sections #use either
#South-North slice at x=1000m
#=

Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=1000.0, orientation=:NS)
plotFlowFieldCrossSection(ControlPlots.plt, wf, A, Zh, Z, vis, 1000.0; orientation=:NS)




loc = 4000.0
Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=loc, orientation=:WE)
plotFlowFieldCrossSection(ControlPlots.plt, wf, A, Zh, Z, vis, loc; orientation=:WE)

=#
#endregion cross-sections



println()
include("calculate_increase_over_baseline.jl")  #to calculate the increase over baseline for the optimised case, compared to a baseline case with no yawing





GC.gc()

#endregion PLOTTING & POSTPROCESSING



println("==========================================")
#notes----------------------------------------------------------------------------------------------------------
#=






















=#


