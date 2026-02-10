#run this file after running test_env/parallel_floridyn.jl
#to see the results of the optimization in a dynamic simulation


x = result.minimizer


    #printing:
    println("-----------")
    println("Optimal yaw change timestamps:")
    println([0.00; x[1:set_num_yaw_changes-1]])
    println()
    println("Optimal yaw angles:")

    for i in 1:set_num_yaw_changes
        println(round.(x[set_num_yaw_changes + (i-1)*wf.nT : set_num_yaw_changes + i*wf.nT - 1] .* 360, digits=1))
    end
    println("-----------")



using Plots

# 1. Setup dimensions
n_intervals = set_num_yaw_changes
n_turbines = wf.nT

# 2. Extract Timestamps
# User logic: x[1:n_intervals-1] are timestamps. 
# We add 0.0 at the start and 1.0 at the end to define the full normalized duration.
t_switches = [0.0; x[1:n_intervals-1]; 1.0]

# 3. Extract and Reshape Yaws
# The yaw data starts at x[n_intervals]
yaw_data_flat = x[n_intervals:end]
# Reshape into (intervals x turbines)
yaws = reshape(yaw_data_flat, n_turbines, n_intervals)' 

# 4. Generate Plot
plt = plot(
    title = "Optimal Turbine Yaw Schedule",
    xlabel = "Normalized Time (0-1)",
    ylabel = "Yaw Angle (Degrees)",
    legend = :outerright,
    grid = true
)

for i in 1:n_turbines
    # We use t_switches for the x-axis. 
    # Because there is one more timestamp than intervals, 
    # we repeat the last yaw value to finish the step plot.
    y_values = yaws[:, i] .* 360.0  # Convert from normalized to degrees
    push!(y_values, y_values[end])  # Match length of t_switches for the step
    
    plot!(plt, t_switches, y_values, 
          label = "Turbine $i", 
          linetype = :steppost, 
          linewidth = 2)
end

display(plt)