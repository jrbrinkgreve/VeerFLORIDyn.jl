#code for plotting the cost function:

#three turbines,plot angles of first 2
#slightly mis-positioned with respect to the wind direction
#can also vary veer magnitude for cost function changes

using Infiltrator


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




#just to get the size
x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)  #start from scratch for baseline
con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)

function construct_yaw_matrix_static!(buffer, x1, x2)
    buffer[:, 2] .= x1
    buffer[:, 3] .= x2 
    buffer[:, 4] .= 180.0

    return nothing
end


function calc_power(md)
    return sum(md.PowerGen) / (wf.nT * sim.n_sim_steps) * 1000.0
end


matrix = zeros(41, 41)
x_array = range(160.0, 200.0, length=41)
y_array = range(160.0, 200.0, length=41)

for (i, x1) in enumerate(x_array)
    for (j, x2) in enumerate(y_array)
        construct_yaw_matrix_static!(con.yaw_data, x1, x2)
        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
        matrix[i, j] = calc_power(md)  # <- whatever your output metric is
    end
end


#wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
#baseline_power_avg = calc_power(md)



