#this code utilises the optimiser in parallel_yaw_optimisation in a receding horizon manner for MPC applications


#=
outline:

if we just make a con.yaw constructor which fixed the past control input, and perturbs the
  wind.dir data, we can easily make an MPC-type controller which has the established flowfield data,
  predicts the best instruction for the future, and then applied the first step!


in pseudocode:


data = load_past_data
weather_prediction = load_weather_prediction
con.yaw_data = optimise(past_info, prediction)
apply_control(con.yaw_data)

repeat

=#

include("MPC_functions.jl")






horizon = 3600 #seconds
#start the time for real-time tracking of the MPC loop
tic()



#write horizon-limited to WindDir.csv, WindVel.csv, WindTI.csv, and set the static floris.veer_gradient later
time = toc()
all_data = load_all_data()
write_horizon_limited_csvs(horizon, time, all_data, path_to_CSVs)




#call initial optimisation with time-selection from 0
#writes a control input
parallel_yaw_optimisation_MPC(init = true)



#"apply" step: visualise it
display_MPC_control()




#when done: call a new optimisation routine with the updated horizon
#needs:
#updated data
#old controller data
#current time 
#





tic()
while time < end_time
    time = get_time()



    all_data = load_all_data()




    write_horizon_limited_csvs(horizon, time, all_data, path_to_CSVs)

    
    
    
    parallel_yaw_optimisation_MPC(init = false)
    
    
    
    
    
    display_MPC_control()



end




















nothing









