using FLORIDyn, Evolutionary, JLD2


#copy this from parallel_floridyn.jl


#optimisation constraints
set_num_yaw_changes = 6 #N
set_max_yaw_rate = 1.0 #deg/s
set_max_yaw_misalignment = 45.0 #deg, for penalising large yaw angles in the cost function, for stability and convergence reasons
#set_num_optimiser_runs = 1  #number of automatic restarts for CMA-ES, unused at the moment
set_sigma0 =  0.02 # 0.01 works well!!          # set to 30% of the search range, and for yaw convergence: first 0.1 for time , then 0.03 for yaws
set_lambda_l1 = 0.0 #1e3  #units: cost PER DEGREE, per turbine, PER SECOND 
                        #relative to the beneficial term average kW per turbine 
                        #typical value 1e3, can play around with this
set_lambda_l1_hard_limit = Inf #a limit on the maximum total yaw change in a simulation, in degrees



@load "workspace.jld2"

include("create_video_files.jl")



