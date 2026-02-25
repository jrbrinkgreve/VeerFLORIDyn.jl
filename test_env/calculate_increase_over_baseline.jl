#this script gives the % increase over baseline performance defined by a wind aligned control strategy



#reset states:
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








x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)  #start from scratch for baseline
con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
baseline_power_avg = sum(md.PowerGen) / (wf.nT * sim.n_sim_steps) * 1000.0




#reset states again:

#reset states:
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

x = result.minimizer  #start from previous result
con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
optimized_power_avg = sum(md.PowerGen) / (wf.nT * sim.n_sim_steps) * 1000.0 #in kW



println("Baseline average power per turbine:  $(round(baseline_power_avg, digits=2)) kW")
println("Optimized average power per turbine: $(round(optimized_power_avg, digits=2)) kW")

increase_over_baseline = (optimized_power_avg - baseline_power_avg) / baseline_power_avg * 100.0

println("Increase over baseline: $(round(increase_over_baseline, digits=2)) %")

