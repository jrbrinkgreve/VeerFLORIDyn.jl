using Timers
using Evolutionary
using OhMyThreads: TaskLocalValue
using Infiltrator
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end


#get the visualisation and settings
_ , vis_file = get_default_project()[2:3]
settings_file = "data/CONTROLTEST_VEER.yaml"   #custom data file with veer specification

vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end

include("../examples/remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

#manually set the initial yaw angles for the optimisation 
con.yaw_data = [sim.start_time:sim.end_time        180 .*  ones(sim.end_time-sim.start_time+1, 9)]

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();


wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);

#run initial conditions
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 

#------------------------------------------------------------------------------------------------------------





function create_fitness(plt, set, wf::WindFarm, wind::Wind, sim, con, vis, floridyn, floris)
    # Create task-local copies of mutable state objects
    # These get deep-copied once per task, ensuring no shared state
    tlv = TaskLocalValue{NamedTuple}() do
        (
            plt = deepcopy(plt),
            set = deepcopy(set),
            wf = deepcopy(wf),
            wind = deepcopy(wind),
            sim = deepcopy(sim),
            floris = deepcopy(floris),
            floridyn = deepcopy(floridyn)
            #need to figure out which input arguments get modified 
            # in-place to determine what needs to be deep-copied
        )
    end
    
    return function cost(x)
        # Get this task's private copies of all mutable state
        con.yaw_data = con.yaw_data = [sim.start_time:sim.end_time        x[1] .*  ones(sim.end_time-sim.start_time+1, 9)]

        state = tlv[]
        # Call runFLORIDyn with task-local copies
        # So each thread modifies its own copies, never shared
        plt_local  = state.plt
        set_local  = state.set
        wf_local = state.wf
        wind_local = state.wind
        sim_local = state.sim
        floris_local = state.floris
        floridyn_local = state.floridyn
        
        # Now safe: each thread has isolated state
        wf, md, mi = run_floridyn(plt_local, set_local, wf_local, wind_local, sim_local, con, vis, floridyn_local, floris_local)
        return -sum(md.PowerGen)
        # Compute fitness from your simulation results
        
    end
end

# Create fitness function (captures all state as closures)
fit_func = create_fitness(plt, set, wf, wind, sim, con, vis, floridyn, floris)

x0 = [185.0]

opts = Evolutionary.Options(
    iterations = 10,
    abstol = 1e-8,
    reltol = 1e-8,
    show_trace = true,
    store_trace = true,
    parallelization = :thread #thread
)

result = Evolutionary.optimize(fit_func, x0, CMAES(), opts)









