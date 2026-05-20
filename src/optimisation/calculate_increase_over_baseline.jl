#this script gives the % increase over baseline performance defined by a wind aligned control strategy
using Plots
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





num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip #skip the first 125 timesteps to avoid startup effects, which is 25 seconds at 5Hz


x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)  #start from scratch for baseline
con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
l1_baseline_yaw_avg = l1_norm_calc(con.yaw_data[:, 2:end]) * (sim.end_time - sim.start_time + 1) 
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
baseline_power_avg = sum(md.PowerGen[(wf.nT*num_timesteps_to_skip+1):end] ) / (wf.nT * (sim.n_sim_steps - num_timesteps_to_skip)) * 1000.0
#feasibility_baseline = cost_func(x0)

# Build time axis (synced to flow sim, skipping initialisation period)
time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step

# Reshape PowerGen into (n_timesteps, nT) and sum across turbines
power_matrix = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)  # (nT, n_sim_steps)
total_power = vec(sum(power_matrix, dims=1))       # kW, one value per timestep
power_turbine_8_trimmed = power_matrix[8, (num_timesteps_to_skip+1):end] * 1000.0  # Power of turbine 8 in kW, trimmed to match time axis


# Slice off the initialisation period
total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]



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



# Extract optimised power timeseries
power_matrix_opt = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
total_power_opt = vec(sum(power_matrix_opt, dims=1))
total_power_opt_trimmed = total_power_opt[(num_timesteps_to_skip+1):end]
power_opt_turbine_8_trimmed = power_matrix_opt[8, (num_timesteps_to_skip+1):end] * 1000.0  # Power of turbine 8 in kW, trimmed to match time axis


plot_P_traces = true

if plot_P_traces
    # Create twin-axis plot: total farm power (left) + turbine 8 power (right)

    plt_power = Plots.plot(
        time_axis, total_power_trimmed,
        xlabel     = "Time (s)",
        ylabel     = "Total Farm Power\n(MW)",
        title      = "Total Farm Power vs Time",
        label      = "Baseline",
        lw         = 2,
        color      = :steelblue,
        legend     = :top,
        legendfontsize = 8,
        grid       = true,
        gridalpha  = 0.4
    )
   Plots.plot!(
        plt_power,
        time_axis, total_power_opt_trimmed,
        label = "Optimised",
        lw    = 2,
        color = :coral
    )

    display(plt_power)
    Plots.savefig(plt_power, "output/power_comparison.pdf")
end


plot_yaw_traces = false




if plot_yaw_traces
    n_turbines = size(con.yaw_data, 2) - 1  # 10 turbines
    time_yaw   = con.yaw_data[num_timesteps_to_skip * sim.time_step + 1:end, 1]
    trimmed_yaw_data = con.yaw_data[num_timesteps_to_skip * sim.time_step + 1:end, 2:end]

    
    # Build a colour palette
    colors = palette(:tab10, n_turbines)
    
    plt_yaw = Plots.plot(
        time_yaw, trimmed_yaw_data[:, 1], 
        xlabel         = "Time (s)",
        ylabel         = "Yaw Angle (°)",
        title          = "Turbine Yaw Angles vs Time",
        label          = "T1",
        lw             = 1.5,
        color          = colors[1],
        legend         = :topright,
        legendfontsize = 7,
        grid           = true,
        gridalpha      = 0.4
    )

    for i in 2:n_turbines
        Plots.plot!(
            plt_yaw,
            time_yaw, trimmed_yaw_data[:, i],
            label = "T$i",
            lw    = 1.5,
            color = colors[i]
        )
    end

    display(plt_yaw)
    Plots.savefig(plt_yaw, "output/yaw_angles.pdf")
end


#=
    # Twin axis for both turbine 8 traces
    ax2 = twinx(plt_power)
    Plots.plot!(
        ax2,
        time_axis, power_turbine_8_trimmed,
        ylabel         = "Turbine 8 Power\n(kW)",
        label          = "Turbine 8 (Baseline)",
        lw             = 2,
        color          = :steelblue,
        linestyle      = :dash,
        legend         = :bottomright,
        legendfontsize = 8,
        grid           = false
    )
    Plots.plot!(
        ax2,
        time_axis, power_opt_turbine_8_trimmed,
        label          = "Turbine 8 (Opt)",
        lw             = 2,
        color          = :coral,
        linestyle      = :dash,
        legendfontsize = 8,
        grid           = false
    )
=#


println("Baseline average power per turbine:  $(round(baseline_power_avg, digits=2)) kW")
println("Optimized average power per turbine: $(round(optimized_power_avg, digits=2)) kW")

energy_increase_over_baseline = (optimized_power_avg - baseline_power_avg) / baseline_power_avg * 100.0

println("Increase over baseline: $(round(energy_increase_over_baseline, digits=2)) %")


println()
println("Baseline average L1 yaw change norm: $(round(l1_baseline_yaw_avg, digits=2)) deg")
println("Optimized average L1 yaw change norm: $(round(l1_optimized_yaw_avg, digits=2)) deg")
println()
#println("Baseline feasibility/direct cost function return: $(round(feasibility_baseline, digits=2))")
#println("Optimized feasibility/direct cost function return: $(round(feasibility_optimized, digits=2))")
#println()




plot_P_percent = false

if plot_P_percent
    P_percentage_increase = (total_power_opt_trimmed - total_power_trimmed) ./ total_power_trimmed * 100.0
    
    # Removed the comma at the end of this line
    P_percent_plot = Plots.plot(
        time_axis, P_percentage_increase,
        xlabel = "Time (s)",
        ylabel = "Increase over Baseline (%)",
        title  = "P% vs Time",
        lw     = 2
    ) 
    
    # Removed the comma at the end of this line
    Plots.plot!(
        P_percent_plot,
        time_axis, ones(length(time_axis)) * energy_increase_over_baseline, ls = :dash, label = "Average Increase"
    )
    
    # Removed the comma at the end of this line
    Plots.plot!(
        P_percent_plot,
        time_axis, zeros(length(time_axis)), ls = :dash, label = "Baseline"
    )

    display(P_percent_plot)
    Plots.savefig(P_percent_plot, "output/power_percentage_increase.pdf")
end



#plotting the increase over baseline:
function plot_optimiser_trace()
    all_points = []
    time_offset = 0.0
    for tr in all_traces
        for record in tr
            t = record.metadata["time"] + time_offset
            v = record.value
            push!(all_points, (t, v))
        end
        time_offset += tr[end].metadata["time"]
    end
    times       = [p[1] for p in all_points]
    vals        = [p[2] for p in all_points]
    running_min = accumulate(min, vals)
    trace_matrix = hcat(times, running_min)

    increase_over_baseline = (-trace_matrix[:, 2] .- baseline_power_avg) ./ baseline_power_avg .* 100.0

    # limit to ±3% range
    ymin, ymax = -1, 2

    plot(trace_matrix[:, 1], increase_over_baseline,
        xlabel        = "Time (s)",
        ylabel        = "Increase over Baseline\n(%)",
        title         = "Optimiser Convergence",
        ylims         = (ymin, ymax),
        yticks        = ymin:0.5:ymax,
        color         = :steelblue,
        linewidth     = 2,
        legend        = false,
        grid          = true,
        gridalpha     = 0.4,
        size          = (700, 350),
        left_margin   = 5Plots.mm,
        right_margin  = 15Plots.mm,
        bottom_margin = 5Plots.mm)
    savefig("output/optimiser_convergence_over_time.pdf")
end

plot_optimiser_trace()




