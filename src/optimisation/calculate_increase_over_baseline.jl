#this script gives the % increase over baseline performance defined by a wind aligned control strategy

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





num_timesteps_to_skip = 125 #skip the first 125 timesteps to avoid startup effects, which is 25 seconds at 5Hz


x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)  #start from scratch for baseline
con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
l1_baseline_yaw_avg = l1_norm_calc(con.yaw_data[:, 2:end]) * (sim.end_time - sim.start_time + 1) 
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
baseline_power_avg = sum(md.PowerGen[(wf.nT*num_timesteps_to_skip+1):end] ) / (wf.nT * (sim.n_sim_steps - num_timesteps_to_skip)) * 1000.0
feasibility_baseline = cost_func(x0)



#reset states again:

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

x = result.minimizer  #start from previous result
con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
l1_optimized_yaw_avg = l1_norm_calc(con.yaw_data[:, 2:end]) * (sim.end_time - sim.start_time + 1) 
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
optimized_power_avg = sum(md.PowerGen[(wf.nT*num_timesteps_to_skip+1):end] ) / (wf.nT * (sim.n_sim_steps - num_timesteps_to_skip)) * 1000.0
feasibility_optimized = cost_func(x)

#plotting the power:



println("Baseline average power per turbine:  $(round(baseline_power_avg, digits=2)) kW")
println("Optimized average power per turbine: $(round(optimized_power_avg, digits=2)) kW")

energy_increase_over_baseline = (optimized_power_avg - baseline_power_avg) / baseline_power_avg * 100.0

println("Increase over baseline: $(round(energy_increase_over_baseline, digits=2)) %")


println()
println("Baseline average L1 yaw change norm: $(round(l1_baseline_yaw_avg, digits=2)) deg")
println("Optimized average L1 yaw change norm: $(round(l1_optimized_yaw_avg, digits=2)) deg")
println()
println("Baseline feasibility/direct cost function return: $(round(feasibility_baseline, digits=2))")
println("Optimized feasibility/direct cost function return: $(round(feasibility_optimized, digits=2))")
println()







