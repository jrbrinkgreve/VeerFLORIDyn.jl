#code for plotting the cost function:

#three turbines,plot angles of first 2
#slightly mis-positioned with respect to the wind direction
#can also vary veer magnitude for cost function changes




include("functions.jl")

#reset states:
# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)


# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();


wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);

#run initial conditions
wf = initSimulation(wf, sim);




#set the situation: 
#=TO DO:
- change turbine location file for 3 in a line or slightly diagonally placed
- double for loop for first and second turbine yaw angles, with third fixed at alignment
- generate fixed control input, and run floridyn for each combination of angles, store the power output in a matrix,
and plot as a contour plot or surface plot
=#

set_num_yaw_changes = 1




x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)  #start from scratch for baseline
con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
baseline_power_avg = sum(md.PowerGen) / (wf.nT * sim.n_sim_steps) * 1000.0



