using Infiltrator



"""
    max_yaw_misalignment(x) -> Float64

Returns the maximum absolute yaw misalignment (degrees) across all turbines
and all timesteps over the full simulation.

Arguments:
  - x               : optimizer decision variable vector (same format as construct_yaw_matrix_dynamic)

The function reconstructs the yaw matrix from x, interpolates wind direction
at every simulation timestep, and computes max |yaw - wind_dir| over all
turbines and timesteps.
"""

function max_yaw_misalignment(x)

    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    # create settings struct with automatic parallel/threading detection
    set = Settings(wind, sim, con, false, false)
    set.enable_veer = true
    set.control_mode = Yaw_Optimisation();
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
    con.yaw_data = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate yaw matrix with zeros to prevent error on 'nothing'
    #run initial conditions
    wf = initSimulation(wf, sim);

    #disable online visualisation
    vis.online = false 
    #construct yaw matrix from x
    yaws = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
    yaws_only = yaws[:, 2:end] #remove time column for misalignment calculation

    #set full wind matrix
    winds = zeros(sim.end_time - sim.start_time + 1)
    fill_wind_dir_buffer!(winds, yaws[:,1], wind.dir)
    

    val = maximum(abs.(yaws_only .- winds))
    arg = argmax(abs.(yaws_only .- winds))
    return val, arg
end












x =3951.13


compare = 3855.04

round((x - compare) / compare * 100, digits=2) 














