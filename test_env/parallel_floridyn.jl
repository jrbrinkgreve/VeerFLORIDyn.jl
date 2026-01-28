using Evolutionary
using OhMyThreads: TaskLocalValue

# to load the sim con plt wf wind floris floridyn set vis objects
include("../examples/veer_main_mini.jl")  



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
        state = tlv[]

        con.yaw_fixed = x[1]
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
        #return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    end
end

# Initialize base objects ONCE before parallelization
plt = plt # your plot object
set = set # your Settings
wf = wf # your WindFarm
wind = wind # your Wind
sim = sim # your simulation state
con = con # your control
vis = vis # your visualization
floridyn = floridyn # your floridyn state
floris = floris # your floris state

# Create fitness function (captures all state as closures)
fit_func = create_fitness(plt, set, wf, wind, sim, con, vis, floridyn, floris)

x0 = [200.0]

opts = Evolutionary.Options(
    iterations = 20,
    abstol = 1e-8,
    reltol = 1e-8,
    show_trace = true,
    store_trace = true,
    parallelization = :thread
)

result = Evolutionary.optimize(fit_func, x0, CMAES(), opts)
