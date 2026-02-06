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
#con.yaw_data = [sim.start_time:sim.end_time        180 .*  ones(sim.end_time-sim.start_time+1, 9)]


# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();


wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);

#run initial conditions
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 

#------------------------------------------------------------------------------------------------------------


function construct_yaw_matrix(x, sim, wf)   #wf to be used later
    yaws =  ones(sim.end_time - sim.start_time + 1)   * x' * 360  #expands into a matrix,
                                                                    # x is in [0,1] range

    #ADD CHECK: yaws must not excees X degrees misalignment from wind to prevent crash


    return   [sim.start_time:sim.end_time    yaws]
end









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
 
        #mapping from x to yaw matrix
        con.yaw_data = construct_yaw_matrix(x, sim, wf)    
        
        
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





#OPTIMISATION PART
#create local copies of vaiables for multithreading


# Create fitness function (captures all state as closures)
fit_func = create_fitness(plt, set, wf, wind, sim, con, vis, floridyn, floris)

x0 = 182/360 * ones(wf.nT)  #182 deg
#x0 = result.minimizer  #start from previous result


#AAAAAAA note for next time: make these limits dependent on the dynamic wind field
lower_bounds = 0.3 * ones(wf.nT)
upper_bounds = 0.7 * ones(wf.nT)


opts = Evolutionary.Options(
    iterations = 10,
    abstol = 1e-8,
    reltol = 1e-8,
    show_trace = true,
    show_every = 1,
    store_trace = true,
    parallelization = :thread  #serial     #thread   
)




#hyperparams
set_lambda_multiplier = 50
set_lambda0 = 2 * round(   (4 + 3 * log(wf.nT)) / 2.0   )   #half the offspring 
set_lambda = Int(set_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))
set_sigma0 = 0.03   # set to 30% of the search range




result =  Evolutionary.optimize(fit_func,
                                BoxConstraints(lower_bounds, upper_bounds),
                                x0, 
                                CMAES(  lambda = set_lambda,
                                        mu = set_mu,
                                        sigma0 = set_sigma0),
                                opts)





#=


#hyperparams
set_lambda_multiplier = 100
set_lambda0 = 2 * round(   (4 + 3 * log(wf.nT)) / 2.0   )
set_lambda = Int(set_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))
set_sigma0 = 0.03   # set to 30% of the search range





result.minimizer = 

no veer:
0.4948785267336715
 0.4976951720463655
 0.5166877606057314
 0.3991364900703976
 0.39732258773717166
 0.42303072369017003
 0.42829439616573045
 0.4354361542932078
 0.44010866886755406






 0.05 veer:
 0.5005065208003683
 0.5062194950628727
 0.4988809864260414
 0.42650958889740204
 0.41829690690115723
 0.4254264869849687
 0.4896144578065447
 0.49324825480570567
 0.4909384118734974




=#




#----------------------------------------------------------------------------------

#PLOTTING

# the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

#create settings struct
#
#set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set = Settings(wind, sim, con, false, false)

set.enable_veer = true
set.control_mode = Yaw_Optimisation();


wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);

#run initial conditions
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 

con.yaw_data = construct_yaw_matrix(result.minimizer, sim, wf)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)










#=
AAAAAAA        
#the idea is now to limit the yaw angles explored by the optimiser: limit them from max.
#+/- 90 deg from current wind direction








I think something is wrong with the multithreading is a bit off,
as the globally best result is not returned while calling get_energy.jl



#next time: IPOP restart strategies for CMA-ES to get out of local minima






=#