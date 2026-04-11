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
settings_file = "data/REALWF_CONTROLTEST_VEER.yaml"   #custom data file with veer specification
vis = Vis(vis_file)
plt=nothing

include("../../examples/remote_plotting.jl")
include("optimisationstructs.jl")
include("functions.jl")

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

#endregion

#------------------------------------------------------------------------------------------------------------------
#OPTIMISATION PART

#sim 
set_num_yaw_changes = 2         #N
set_max_yaw_misalignment = 45.0 #deg, for penalising large yaw angles in the cost function, for stability and convergence reasons
set_lambda_l1 = 0            #1e3  #units: cost PER DEGREE, per turbine, PER SECOND #relative to the beneficial term average kW per turbine #typical value 1e3, can play around with this
set_lambda_l1_hard_limit = Inf  #a limit on the maximum total yaw change in a simulation, in degrees
set_max_yaw_rate = 1.0          #deg/s
set_objective = totalEnergyObjective      #totalEnergyObjective or powerTrackingObjective 

#optimiser convergence/ fidelity  
set_cmaes_lambda_multiplier = 4    #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
set_num_optimiser_runs = 3        #number of automatic restarts for CMA-ES
set_iterations = 50                #number of iterations for CMAES
set_sigma0 = 0.05                  # 0.05 for time optim, 0.01 works well in second run for yaws!!     # set to 30% of the search range, and for yaw convergence: first 0.1 for time , then 0.03 for yaws
set_sigma0_secondary = 0.01         #for second run with yaws, to reduce the search area and converge faster, can play around with this


#initial state
x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)   
#x0 = result.minimizer           #start from previous result also possible

#region fold CMAES prep


#efficient struct passing
opt_set = OptimisationSettings(
    set_num_yaw_changes, set_num_optimiser_runs, set_max_yaw_rate, set_max_yaw_misalignment, set_lambda_l1, set_lambda_l1_hard_limit, set_objective
)

#hyperparams
#set_lambda0 = 10
set_lambda0 = 2 * round(   (4 + 3 * log(wf.nT * (set_num_yaw_changes-1)))      / 2.0   )   #half the offspring 
set_lambda = Int(set_cmaes_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))



# set cost function
cost_func = parallel_costfunction(plt, set, wf, wind, sim, con, vis, floridyn, floris, opt_set)


#opts
opts = Evolutionary.Options( 
    iterations = set_iterations,
    abstol = 1e-8,
    reltol = 1e-8,          #1e-8
    show_trace = true,
    show_every = 5, 
    store_trace = true,
    parallelization = :thread  #serial   #thread           #multithreading hehe
)



#constraints
#limit normalised time to [0,1], yaws are free after adding try/catch
lower_bounds = zeros(set_num_yaw_changes-1)
upper_bounds = ones(set_num_yaw_changes-1)
#endregion CMAES prep


println()
println("Starting CMAES run 1 / $set_num_optimiser_runs")
println("σ_0 = $set_sigma0, λ = $set_lambda, μ = $set_mu")
println()
all_traces = [] #to store traces from multiple runs 

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


#secondary runs with other options as well
for run in 2:set_num_optimiser_runs

    #hyperparams
    local_set_lambda_multiplier = set_cmaes_lambda_multiplier #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
    local_set_lambda0 = set_lambda0 
    local_set_lambda = set_lambda
    local_set_mu = set_mu

    #reduce sigma to stay in local region
    local_set_sigma0 = set_sigma0_secondary

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
end




end_time = time()
total_time = end_time - begin_time
println()
println("Total optimization time: $(round(total_time, digits=2)) seconds")

#endregion CMAES

#------------------------------------------------------------------------------------------------------------------

#region PLOTTING & POSTPROCESSING
#=
#convergence plot
threshold = 1e5
p = Plots.plot(title="CMA-ES Convergence (all runs)", xlabel="Iteration", ylabel="Cost")

let offset = 0
    for (i, trace) in enumerate(all_traces)
        filtered = [(t.iteration, t.value) for t in trace if t.value < threshold]
        iters    = [x[1] + offset for x in filtered]
        values   = [x[2]          for x in filtered]
        Plots.plot!(p, iters, values, lw=2, label="Run $i")
        offset = isempty(iters) ? offset : iters[end]
    end
end
display(p)
=#


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


#plot flowfield
#construct_yaw_matrix_dynamic!(con.yaw_data, result.minimizer, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
construct_yaw_matrix_dynamic!(con.yaw_data, result.minimizer, sim, wf, opt_set)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

#top-view
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)







#region cross-sections #use either
#South-North slice at x=1000m
#=

Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=1000.0, orientation=:NS)
plotFlowFieldCrossSection(ControlPlots.plt, wf, A, Zh, Z, vis, 1000.0; orientation=:NS)

=#

#West-East slice at y=1500m
#=

Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=2000.0, orientation=:WE)
plotFlowFieldCrossSection(ControlPlots.plt, wf, A, Zh, Z, vis, 2000.0; orientation=:WE)

=#
#endregion cross-sections



println()
include("calculate_increase_over_baseline.jl")  #to calculate the increase over baseline for the optimised case, compared to a baseline case with no yawing



GC.gc()

#endregion PLOTTING & POSTPROCESSING


#notes----------------------------------------------------------------------------------------------------------
#=






















=#






