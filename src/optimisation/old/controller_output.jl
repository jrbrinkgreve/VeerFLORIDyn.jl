#run this file after running test_env/parallel_floridyn.jl
#to see the results of the optimization in a dynamic simulation


yaws = construct_yaw_matrix_dynamic(result.minimizer, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
yaws = yaws[:, 2:end]  #remove time column

using Plots

# Plot the yaw angles over time for each turbine
plt = plot(
    title = "Optimal Turbine Yaw Schedule",
    xlabel = "Time (s)",
    ylabel = "Yaw Angle (Degrees)",
    legend = :outerright,
    grid = true
)

for i in 1:wf.nT
    plot!(plt, yaws[:, i], label = "Turbine $i", linewidth = 2)
end

display(plt)
