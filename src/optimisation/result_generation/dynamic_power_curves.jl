using Statistics
using Plots
using JLD2
using Plots
using Interpolations
using LaTeXStrings
include("../functions.jl")



color_palette = [:coral, :mediumseagreen, :blue, :red, :purple]


function lambda_curves()

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

# Build time axis (synced to flow sim, skipping initialisation period)
time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step

# Reshape PowerGen into (n_timesteps, nT) and sum across turbines
power_matrix = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)  # (nT, n_sim_steps)
total_power = vec(sum(power_matrix, dims=1))       # kW, one value per timestep
power_turbine_8_trimmed = power_matrix[8, (num_timesteps_to_skip+1):end] * 1000.0

# Slice off the initialisation period
total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]

# Bold colour palette
lambda_colors = color_palette

# 1. Initialize the plot with the Baseline trace
plt_power = Plots.plot(
        time_axis, ones(length(time_axis)),
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power \n(-)",
        #title          = "Baseline-normalised Farm Power vs Time: λ-regularisation",
        label          = "Baseline",
        lw             = 3,
        color          = :black,
        legend         = :outerright,
        legendfontsize = 15,
        tickfontsize   = 12,
        guidefontsize  = 12,
        grid           = true,
        gridalpha      = 0.4
    )

# 2. Define the lambda values you want to test
lambdas = [0.0, 2000.0, 7000.0, 20000.0]

# 3. Loop over each lambda, load the result, run the simulation, and append to the plot
for (i, λ) in enumerate(lambdas)
    println("Running simulation for lambda = $λ...")
    
    # Load the specific minimizer for this lambda
    path = "src/optimisation/optim_data/minimisers/lambda_reg_$λ.jld2"
    result = load(path)["result"]
    x = result.minimizer

    # Reset states completely for a fresh simulation
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    # Create settings struct 
    set = Settings(wind, sim, con, false, false)
    set.enable_veer  = true
    set.control_mode = Yaw_Optimisation()

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

    # Run initial conditions
    wf = initSimulation(wf, sim)

    # Apply the optimized yaw layout from JLD2
    con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
    
    l1_optimized_yaw_avg = l1_norm_calc(con.yaw_data[:, 2:end]) * (sim.end_time - sim.start_time + 1) 
    
    # Run the actual simulation
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

    # Extract optimised power timeseries
    power_matrix_opt = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
    total_power_opt = vec(sum(power_matrix_opt, dims=1))
    total_power_opt_trimmed = total_power_opt[(num_timesteps_to_skip+1):end]
    total_power_opt_trimmed_normalised = total_power_opt_trimmed ./ total_power_trimmed 

    Plots.plot!(
        plt_power,
        time_axis, total_power_opt_trimmed_normalised,
        label = "(λ = $λ)",
        lw    = 3,
        color = lambda_colors[i],
        alpha = 1
    )

end


Plots.plot!(
        plt_power,
        size          = (1250, 350),
        left_margin   = 10Plots.mm,
        right_margin  = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.8,
        tickfontsize  = 12,
        guidefontsize = 12,
        legendfontsize = 15,
        ylims  = (0.85, 1.15),
        yticks = 0.85:0.05:1.15,
        xticks = 500:500:4000
    )

display(plt_power)
Plots.savefig(plt_power, "output/power_comparison_lambdas.pdf")
println("Plot complete. Saved to output/power_comparison_lambdas.pdf")


end







function sparsity_curves()

    # Reset states
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    set = Settings(wind, sim, con, false, false)
    set.enable_veer  = true
    set.control_mode = Yaw_Optimisation()

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    wf = initSimulation(wf, sim)

    num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip

    x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
    con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

    time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step

    power_matrix      = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
    total_power       = vec(sum(power_matrix, dims=1))
    total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]

    # Bold colour palette
    sparsity_colors = color_palette

    plt_sparsity = Plots.plot(
        time_axis, ones(length(time_axis)),
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power\n(-)",
        #title          = "Baseline-normalised Farm Power vs Time: Sparsity constraints",
        label          = "Baseline",
        lw             = 3,
        color          = :black,
        legend         = :outerright,
        tickfontsize   = 12,
        guidefontsize  = 12,
        legendfontsize = 15,
        grid           = true,
        gridalpha      = 0.4
    )

    sparsities = [2, 3, 4, 6]

    for (i, s) in enumerate(sparsities)
        println("Running simulation for sparsity = $s...")
        set_num_yaw_changes = s

        path   = "src/optimisation/optim_data/minimisers/sparsity_$s.jld2"
        result = load(path)["result"]
        x      = result.minimizer
        
        wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

        set = Settings(wind, sim, con, false, false)
        set.enable_veer  = true
        set.control_mode = Yaw_Optimisation()

        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        wf = initSimulation(wf, sim)

        con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
        
        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

        power_matrix_opt            = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
        total_power_opt             = vec(sum(power_matrix_opt, dims=1))
        total_power_opt_trimmed     = total_power_opt[(num_timesteps_to_skip+1):end]
        total_power_opt_normalised  = total_power_opt_trimmed ./ total_power_trimmed

        Plots.plot!(
            plt_sparsity,
            time_axis, total_power_opt_normalised,
            label         = "(s = $s)",
            lw            = 3,
            color         = sparsity_colors[i],
            alpha         = 1,
            tickfontsize  = 12,
            guidefontsize = 12,
            legendfontsize = 15,
        )
    end

    Plots.plot!(
        plt_sparsity,
        size          = (1250, 350),
        left_margin   = 10Plots.mm,
        right_margin  = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.8,
        ylims         = (0.85, 1.15),
        yticks        = 0.85:0.05:1.15,
        xticks        = 500:500:4000
    )

    display(plt_sparsity)
    Plots.savefig(plt_sparsity, "output/power_comparison_sparsity.pdf")
    println("Plot complete. Saved to output/power_comparison_sparsity.pdf")

end



function yaw_limit_curves()

    # Reset states
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    set = Settings(wind, sim, con, false, false)
    set.enable_veer  = true
    set.control_mode = Yaw_Optimisation()

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    wf = initSimulation(wf, sim)

    num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip

    x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
    con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

    time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step

    power_matrix        = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
    total_power         = vec(sum(power_matrix, dims=1))
    total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]

    # Bold colour palette
    yaw_colors = color_palette

    plt_yaw_lim = Plots.plot(
        time_axis, ones(length(time_axis)),
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power\n(-)",
        #title          = "Baseline-normalised Farm Power vs Time: Yaw misalignment limits",
        label          = "Baseline",
        lw             = 3,
        color          = :black,
        legend         = :outerright,
        legendfontsize = 15,
        tickfontsize   = 12,
        guidefontsize  = 12,
        grid           = true,
        gridalpha      = 0.8
    )

    yaw_limits = [10.0, 15.0, 20.0, 25.0, Inf]

    for (i, yaw_lim) in enumerate(yaw_limits)
        println("Running simulation for yaw limit = $yaw_lim...")

        path   = "src/optimisation/optim_data/minimisers/yaw_lim_$yaw_lim.jld2"
        result = load(path)["result"]
        x      = result.minimizer

        wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

        set = Settings(wind, sim, con, false, false)
        set.enable_veer  = true
        set.control_mode = Yaw_Optimisation()

        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        wf = initSimulation(wf, sim)

        con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)

        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

        power_matrix_opt           = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
        total_power_opt            = vec(sum(power_matrix_opt, dims=1))
        total_power_opt_trimmed    = total_power_opt[(num_timesteps_to_skip+1):end]
        total_power_opt_normalised = total_power_opt_trimmed ./ total_power_trimmed

        Plots.plot!(
            plt_yaw_lim,
            time_axis, total_power_opt_normalised,
            label         = "Limit = $(yaw_lim)°",
            lw            = 3,
            color         = yaw_colors[i],
            alpha         = 1,
            tickfontsize  = 12,
            guidefontsize = 12,
            legendfontsize = 15,
        )
    end

    Plots.plot!(
        plt_yaw_lim,
        size          = (1250, 350),
        left_margin   = 10Plots.mm,
        right_margin  = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        top_margin    = 10Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.8,
        ylims         = (0.85, 1.15),
        yticks        = 0.85:0.05:1.15,
        xticks        = 500:500:4000
    )

    display(plt_yaw_lim)
    Plots.savefig(plt_yaw_lim, "output/power_comparison_yaw_limits.pdf")
    println("Plot complete. Saved to output/power_comparison_yaw_limits.pdf")

end


function veer_curves()

    # Only need setup here to compute time_axis — no baseline power run needed
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    set = Settings(wind, sim, con, false, false)
    set.enable_veer  = true
    set.control_mode = Yaw_Optimisation()

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    wf = initSimulation(wf, sim)

    num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip
    time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step

    # Bold colour palette
    veer_colors = color_palette

    plt_veer_grad = Plots.plot(
        time_axis, ones(length(time_axis)),
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power\n(-)",
        #title          = "Baseline-normalised Farm Power vs Time: Veer gradients",
        label          = "All Baselines",
        lw             = 3,
        color          = :black,
        legend         = :outerright,
        legendfontsize = 15,
        tickfontsize   = 12,
        guidefontsize  = 12,
        grid           = true,
        gridalpha      = 0.4
    )

    veer_values = [0, 0.05, 0.1, 0.15, 0.2]

    for (i, veer_grad) in enumerate(veer_values)
        println("Running simulation for veer gradient = $veer_grad...")

        # --- Baseline: x0 (no optimisation), WITH matched veer gradient ---
        wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

        set = Settings(wind, sim, con, false, false)
        set.enable_veer  = true
        set.control_mode = Yaw_Optimisation()

        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        floris.veer_gradient = veer_grad          # ← matched veer for baseline
        wf = initSimulation(wf, sim)

        x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
        con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

        time_axis           = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step
        power_matrix        = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
        total_power         = vec(sum(power_matrix, dims=1))
        total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]

        # --- Optimised: load minimiser, WITH matched veer gradient ---
        path   = "src/optimisation/optim_data/minimisers/veer_$veer_grad.jld2"
        result = load(path)["result"]
        x      = result.minimizer

        wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

        set = Settings(wind, sim, con, false, false)
        set.enable_veer  = true
        set.control_mode = Yaw_Optimisation()

        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        floris.veer_gradient = veer_grad          # ← matched veer for optimised run
        wf = initSimulation(wf, sim)

        con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
        wf, md, mi   = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

        power_matrix_opt           = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
        total_power_opt            = vec(sum(power_matrix_opt, dims=1))
        total_power_opt_trimmed    = total_power_opt[(num_timesteps_to_skip+1):end]
        total_power_opt_normalised = total_power_opt_trimmed ./ total_power_trimmed  # both at same veer level

        Plots.plot!(
            plt_veer_grad,
            time_axis, total_power_opt_normalised,
            label          = "$(veer_grad)°/m veer",
            lw             = 3,
            color          = veer_colors[i],
            alpha          = 1,
            tickfontsize   = 12,
            guidefontsize  = 12,
            legendfontsize = 12,
        )
    end

    Plots.plot!(
        plt_veer_grad,
        size          = (1250, 350),
        left_margin   = 10Plots.mm,
        right_margin  = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        top_margin    = 5Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.4,
        ylims         = (0.95, 1.25),
        yticks        = 0.85:0.05:1.30,
        xticks        = 500:500:4000
    )

    display(plt_veer_grad)
    Plots.savefig(plt_veer_grad, "output/power_comparison_veer_gradients.pdf")
    println("Plot complete. Saved to output/power_comparison_veer_gradients.pdf")

end


function veer_no_veer_setpoint_comparisons()

    control_setpoint_comparison_plot = Plots.plot()

    # Bold colour palette
    setpoint_colors = color_palette

    veer_values = [0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2]

    for (i, veer_grad) in enumerate(veer_values)

        # Power values for no-veer setpoints in veered cases
        wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

        path   = "src/optimisation/optim_data/minimisers/veer_0.0.jld2"
        result = load(path)["result"]
        x      = result.minimizer

        set = Settings(wind, sim, con, false, false)
        set.enable_veer  = true
        set.control_mode = Yaw_Optimisation()

        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        wf = initSimulation(wf, sim)
        floris.veer_gradient = veer_grad

        num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip

        x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
        con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

        time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step
    
        power_matrix        = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
        total_power         = vec(sum(power_matrix, dims=1))
        total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]
        a = mean(total_power_trimmed)

        # Load proper control instruction
        println("Running simulation for veer gradient = $veer_grad...")

        path   = "src/optimisation/optim_data/minimisers/veer_$veer_grad.jld2"
        result = load(path)["result"]
        x      = result.minimizer

        wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

        set = Settings(wind, sim, con, false, false)
        set.enable_veer  = true
        set.control_mode = Yaw_Optimisation()

        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        floris.veer_gradient = veer_grad
        wf = initSimulation(wf, sim)
  
        con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)

        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

        power_matrix_opt           = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
        total_power_opt            = vec(sum(power_matrix_opt, dims=1))
        total_power_opt_trimmed    = total_power_opt[(num_timesteps_to_skip+1):end]
        total_power_opt_normalised = total_power_opt_trimmed ./ total_power_trimmed

        # Quantifying losses
        println("Power increase percentage:")
        increase_percentage = (sum(total_power_opt_trimmed) - sum(total_power_trimmed)) / sum(total_power_trimmed) * 100
        println(a * 100, "kW")
        println()

        Plots.plot!(
            control_setpoint_comparison_plot,
            time_axis, total_power_opt_normalised,
            label         = "$(veer_grad)°/m veer",
            lw            = 3,
            color         = setpoint_colors[i],
            alpha         = 1,
            tickfontsize  = 12,
            guidefontsize = 12,
            legendfontsize = 12,
            legend         = :outerright,
            grid           = true,
            gridalpha      = 0.4,
            xlabel         = "Time (s)",
            ylabel         = "Normalised Power\n(-)",
            title          = "Veer-Unaware Wake Steering Controller in Veered Conditions:"
        )

    end
    Plots.plot!(
        control_setpoint_comparison_plot,
        size          = (1250, 350),
        left_margin   = 10Plots.mm,
        right_margin  = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.4,
        ylims         = (0.90, 1.10),
        yticks        = 0.85:0.05:1.15,
        xticks        = 500:500:4000
    )

    display(control_setpoint_comparison_plot)
    Plots.savefig(control_setpoint_comparison_plot, "output/control_setpoints_comparison.pdf")
    println("Plot complete. Saved to output/control_setpoints_comparison.pdf")

end



function plot_first_turbine_yaws()
    veer_values = [0.0, 0.05, 0.1, 0.15, 0.2]

    # Bold colour palette
    yaw_colors = [:red, :blue, :green, :darkorange, :purple]

    turbine_yaws_plot = Plots.plot()
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    set = Settings(wind, sim, con, false, false)
    set.enable_veer  = true
    set.control_mode = Yaw_Optimisation()

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
 
    wf = initSimulation(wf, sim)

    num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip

    x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
    con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

    for (i, veer_grad) in enumerate(veer_values)
        path   = "src/optimisation/optim_data/minimisers/veer_$(veer_grad).jld2"
        result = load(path)["result"]
        x      = result.minimizer
        yaws   = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)

        wind_dir_interp = interpolate((wind.dir[:, 1],), wind.dir[:, 2], Gridded(Linear()))
        wind_dir_at_yaw = wind_dir_interp.(yaws[:, 1])

        yaw_misalignment = yaws[:, 2] .- wind_dir_at_yaw
        Plots.plot!(
            turbine_yaws_plot,
            yaws[:, 1], yaw_misalignment,
            label          = "$(veer_grad)°/m veer",
            lw             = 3,
            color          = yaw_colors[i],
            alpha          = 1,
            tickfontsize   = 12,
            guidefontsize  = 12,
            legendfontsize = 12,
            legend         = :outerright,
            grid           = true,
            gridalpha      = 0.4,
            xlabel         = "Time (s)",
            ylabel         = "Yaw misalignment (°)",
            title          = "Dynamic Yaw Misalignment of T1 across veer levels"
        )
    end

    Plots.plot!(
        turbine_yaws_plot,
        size          = (1250, 350),
        left_margin   = 10Plots.mm,
        right_margin  = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.4,
        ylims         = (-30, 30),
        xlims         = (500, 4100),
        yticks        = -30:5:30,
        xticks        = 500:500:4000,
    )

    display(turbine_yaws_plot)
    Plots.savefig(turbine_yaws_plot, "output/veered_first_turbine_yaws.pdf")
    println("Plot complete. Saved to output/veered_first_turbine_yaws.pdf")
end



function plot_first_turbine_yaws_absolute()
    veer_values = [0.0, 0.05, 0.1, 0.15, 0.2]

    # Bold colour palette
    abs_yaw_colors = color_palette

    turbine_yaws_plot = Plots.plot()
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)
    set = Settings(wind, sim, con, false, false)
    set.enable_veer  = true
    set.control_mode = Yaw_Optimisation()
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    wf = initSimulation(wf, sim)
    num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip
    x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
    con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

    for (i, veer_grad) in enumerate(veer_values)
        path   = "src/optimisation/optim_data/minimisers/veer_$(veer_grad).jld2"
        result = load(path)["result"]
        x      = result.minimizer
        yaws   = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)

        wind_dir_interp  = interpolate((wind.dir[:, 1],), wind.dir[:, 2], Gridded(Linear()))
        wind_dir_at_yaw  = wind_dir_interp.(yaws[:, 1])

        if i == 1
            Plots.plot!(
                turbine_yaws_plot,
                yaws[:, 1], wind_dir_at_yaw,
                label     = "Wind Direction",
                color     = :black,
                lw        = 3,
                alpha     = 1.0,
                linestyle = :dash,
            )
        end

        Plots.plot!(
            turbine_yaws_plot,
            yaws[:, 1], yaws[:, 2],
            label          = "veer $(veer_grad)°/m",
            color          = abs_yaw_colors[i],
            lw             = 3,
            alpha          = 1.0,
            tickfontsize   = 12,
            guidefontsize  = 12,
            legendfontsize = 11,
            legend         = :outerright,
            grid           = true,
            gridalpha      = 0.4,
            xlabel         = "Time (s)",
            ylabel         = "Absolute angle (°)",
            title          = "Dynamic Yaw & Wind Direction of T1 across veer levels"
        )
    end

    Plots.plot!(
        turbine_yaws_plot,
        size          = (1250, 350),
        left_margin   = 10Plots.mm,
        right_margin  = 10Plots.mm,
        bottom_margin = 10Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.4,
        ylims         = (170, 230),
        xlims         = (500, 4100),
        xticks        = 500:500:4000,
    )

    display(turbine_yaws_plot)
    Plots.savefig(turbine_yaws_plot, "output/veered_first_turbine_yaws_absolute.pdf")
    println("Plot complete. Saved to output/veered_first_turbine_yaws_absolute.pdf")
end



#lambda_curves()
#sparsity_curves()
#yaw_limit_curves()
veer_curves()
#veer_no_veer_setpoint_comparisons()
#plot_first_turbine_yaws()  
#plot_first_turbine_yaws_absolute()