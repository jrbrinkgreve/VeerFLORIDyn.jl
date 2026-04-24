using FLORIDyn, Evolutionary, JLD2
using ControlPlots
using DistributedNext

#copy this from parallel_floridyn.jl


include("../../examples/remote_plotting.jl")
include("optimisationstructs.jl")
include("functions.jl")



#sim 
set_num_yaw_changes = 4         #N
set_max_yaw_misalignment = 30.0 #deg, for penalising large yaw angles in the cost function, for stability and convergence reasons
set_lambda_l1 = 0.0            #1e3  #units: cost PER DEGREE, per turbine, PER SECOND #relative to the beneficial term average kW per turbine #typical value 1e3, can play around with this
set_lambda_l1_hard_limit = Inf  #a limit on the maximum total yaw change in a simulation, in degrees
set_max_yaw_rate = 1.0          #deg/s
set_objective = totalEnergyObjective   #totalEnergyObjective or powerTrackingObjective al
set_num_timesteps_to_skip = 125     #skip the first N timesteps for wake effects to propagate, approx time between wake interactions



@load "src/optimisation/optim_data/workspace.jld2"

include("create_video_files.jl")



