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
 

        con.yaw_data = construct_yaw_matrix(x, sim, wf)

        #mapping from x to yaw matrix
        #con.yaw_data = construct_yaw_matrix(x, sim, wf)    
        
        
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




"""
 function construct_yaw_matrix()
        #construct yaw matrix for optimisation, reducing the number of free variables for the optimiser


"""

function construct_yaw_matrix(x, sim, wf)
    yaws =  ones(sim.end_time - sim.start_time + 1)   * x' * 360  #expands into a matrix,
                                                                    # x is in [0,1] range

    return   [sim.start_time:sim.end_time    yaws]
end





#OPTIMISATION PART
#create local copies of vaiables for multithreading


# Create fitness function (captures all state as closures)
fit_func = create_fitness(plt, set, wf, wind, sim, con, vis, floridyn, floris)





x0 = 0.5 * 360 * ones(wf.nT)

# ============================================================
# IPOP-CMA-ES Setup
# ============================================================

set_lambda_multiplier = 5.0
set_lambda0 = 2.0 * round((4 + 3 * log(wf.nT)) / 2.0)
set_lambda = Int(set_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))
set_sigma0 = 0.25


# 
    # Options for this restart
    opts = Evolutionary.Options(
        iterations=iterations_per_restart,
        abstol=1e-8,
        reltol=1e-8,
        show_trace=true,
        store_trace=true,
        parallelization=:thread
    )
    
  
    result = Evolutionary.optimize(
        fit_func,
        x0,
        CMAES(
            lambda=set_lambda,
            mu=set_mu,
            sigma0=set_sigma0
        ),
        opts
    )
    result



#plotting the final result:



