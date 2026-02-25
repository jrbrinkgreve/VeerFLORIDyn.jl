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
plt=nothing
#=
vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end
=#

include("../examples/remote_plotting.jl")
include("functions.jl")

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
#OPTIMISATION PART


#number of yaw changes allowed
set_num_yaw_changes = 4 #N
set_max_yaw_rate = 1.0 #deg/s
set_max_yaw_misalignment = 85.0 #deg, for penalising large yaw angles in the cost function, for stability and convergence reasons
#set_num_optimiser_runs = 1  #number of automatic restarts for CMA-ES, unused at the moment

# set cost function
cost_func = parallel_costfunction(plt, set, wf, wind, sim, con, vis, floridyn, floris)




#initial state
#x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
x0 = result.minimizer  #start from previous result





#hyperparams
set_lambda_multiplier = 3 #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
set_lambda0 = 2 * round(   (4 + 3 * log(wf.nT * (set_num_yaw_changes-1)))      / 2.0   )   #half the offspring 
set_lambda = Int(set_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))
set_sigma0 =  0.01  # set to 30% of the search range, and for yaw convergence: first 0.1 for time , then 0.03 for yaws



#opts
opts = Evolutionary.Options( 
    iterations = 100,
    abstol = 1e-8,
    reltol = 1e-8,
    show_trace = true,
    show_every = 1, 
    store_trace = true,
    parallelization = :thread  #serial   #thread     #multithreading hehe
)


#constraints
#limit normalised time to [0,1], yaws are free after adding try/catch
lower_bounds = vcat(zeros(set_num_yaw_changes-1))   #,    -10 * ones(set_num_yaw_changes * wf.nT))
upper_bounds = vcat(ones(set_num_yaw_changes-1))  #,      10 * ones(set_num_yaw_changes * wf.nT))




@time result =  Evolutionary.optimize(cost_func,
                                BoxConstraints(lower_bounds, upper_bounds),
                                x0, 
                                CMAES(  lambda = set_lambda,
                                        mu = set_mu,
                                        sigma0 = set_sigma0),
                                opts)


















#----------------------------------------------------------------------------------

#PLOTTING



# the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)


#set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();

#run initial conditions
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 


#plot flowfield
con.yaw_data = construct_yaw_matrix_dynamic(result.minimizer, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)


#include("get_energy.jl")  #for energy tests
#include("controller_output.jl")   #to see and plot controller angles






#=
AAAAAAA        
#the idea is now to limit the yaw angles explored by the optimiser: limit them from max.
#+/- 90 deg from current wind direction








I think something is wrong with the multithreading is a bit off,
as the globally best result is not returned while calling get_energy.jl

#include automatic IPOP restart strategies for CMA-ES to get out of local minima?

#----------------------------------------------------------------------------------

also some more notes:
20 feb 2026
individual turbine switching has a lower minimum but finding it is too hard for the optimiser,
so we will stick with simultaneous switching for now, and can report about this in the thesis



next steps are minimizing the l1 norm of the yaw actuation on top of the power maximisation. 
the code provides an interface for this via the cost function:
    max_yaw_misalignment is a way of penalising large yaw angles (also for stability)
    set_


=#