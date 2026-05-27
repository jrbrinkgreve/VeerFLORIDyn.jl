using Plots
using JLD2
using Plots
include("../functions.jl")



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
#feasibility_baseline = cost_func(x0)

# Build time axis (synced to flow sim, skipping initialisation period)
time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step

# Reshape PowerGen into (n_timesteps, nT) and sum across turbines
power_matrix = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)  # (nT, n_sim_steps)
total_power = vec(sum(power_matrix, dims=1))       # kW, one value per timestep
power_turbine_8_trimmed = power_matrix[8, (num_timesteps_to_skip+1):end] * 1000.0  # Power of turbine 8 in kW, trimmed to match time axis


# Slice off the initialisation period
total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]


# Assuming you already have `time_axis` and `total_power_trimmed` (Baseline) 
# computed and available in your workspace.

# 1. Initialize the plot with the Baseline trace
plt_power = Plots.plot(
        time_axis, ones(length(time_axis)),  # Start with a dummy trace (will be overwritten)
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power \n(-)",
        title          = "Baseline-normalised Farm Power vs Time",
        label          = "Baseline",
        lw             = 2,
        color          = :black,
        legend         = :bottomright,
        legendfontsize = 10,
        tickfontsize = 12,
        guidefontsize = 12,
        grid           = true,
        gridalpha      = 0.4
    )

# 2. Define the lambda values you want to test
lambdas = [0.0, 2000.0, 7000.0,  20000.0]

# 3. Loop over each lambda, load the result, run the simulation, and append to the plot
for λ in lambdas
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
    
    # (Optional) calculate L1 norm if you need it logged/printed
    l1_optimized_yaw_avg = l1_norm_calc(con.yaw_data[:, 2:end]) * (sim.end_time - sim.start_time + 1) 
    
    # Run the actual simulation
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

    # Extract optimised power timeseries
    power_matrix_opt = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
    total_power_opt = vec(sum(power_matrix_opt, dims=1))
    total_power_opt_trimmed = total_power_opt[(num_timesteps_to_skip+1):end]
    total_power_opt_trimmed_normalised = total_power_opt_trimmed ./ total_power_trimmed 
    
    
    
    
    # Calculate avg power if needed for logging
    # optimized_power_avg = sum(total_power_opt_trimmed) / (wf.nT * (sim.n_sim_steps - num_timesteps_to_skip)) * 1000.0

    # Append the new trace to the existing plot using Plots.plot! (notice the !)
      Plots.plot!(
        plt_power,
        time_axis, total_power_opt_trimmed_normalised,
        label = "(λ = $λ)",
        lw    = 2,
        alpha = 0.8
    )

end


Plots.plot!(
        plt_power,
        size          = (1000, 350),
        left_margin   = 7Plots.mm,
        right_margin  = 15Plots.mm,
        bottom_margin = 5Plots.mm,
        xgrid         = true,
        ygrid         = true,
        gridalpha     = 0.8,
        tickfontsize=12,
        guidefontsize=12,
        legendfontsize=10,
        ylims = (0.85, 1.15),
        yticks = 0.85:0.05:1.15,
        xticks = 500:500:4000

    )


# 4. Display and save the final composed plot
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

    plt_sparsity = Plots.plot(
        time_axis, ones(length(time_axis)),
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power\n(-)",
        title          = "Baseline-normalised Farm Power vs Time",
        label          = "Baseline",
        lw             = 2,
        color          = :black,
        legend         = :bottomright,
        tickfontsize=12,
        guidefontsize=12,
        legendfontsize=10,
        grid           = true,
        gridalpha      = 0.4
    )

    sparsities = [2, 3, 4, 5, 6] # fill in your sparsity values here

    for s in sparsities
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
            label = "(s = $s)",
            lw    = 2,
            alpha = 0.8,
            tickfontsize=12,
            guidefontsize=12,
            legendfontsize=10,
        )
    end

    Plots.plot!(
        plt_sparsity,
        size          = (1000, 350),
        left_margin   = 7Plots.mm,
        right_margin  = 15Plots.mm,
        bottom_margin = 5Plots.mm,
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

    plt_yaw_lim = Plots.plot(
        time_axis, ones(length(time_axis)),
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power\n(-)",
        title          = "Baseline-normalised Farm Power vs Time",
        label          = "Baseline",
        lw             = 2,
        color          = :black,
        legend         = :bottomright,
        legendfontsize = 10,
        tickfontsize=12,
        guidefontsize=12,
        grid           = true,
        gridalpha      = 0.8
    )

    yaw_limits = [10.0 15.0 20.0 25.0 Inf] # fill in your yaw limit values here

    for yaw_lim in yaw_limits
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
            label = "(Yaw limit = $(yaw_lim)°)",
            lw    = 2,
            alpha = 0.8,
            tickfontsize=12,
            guidefontsize=12,
            legendfontsize=10,
        )
    end

    Plots.plot!(
        plt_yaw_lim,
        size          = (1000, 350),
        left_margin   = 7Plots.mm,
        right_margin  = 15Plots.mm,
        bottom_margin = 5Plots.mm,
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

    # Reset states
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    set = Settings(wind, sim, con, false, false)
    set.enable_veer  = true
    set.control_mode = Yaw_Optimisation()

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    floris.veer_gradient =
    wf = initSimulation(wf, sim)

    num_timesteps_to_skip = opt_set.set_num_timesteps_to_skip

    x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
    con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

    time_axis = sim.start_time .+ (num_timesteps_to_skip:sim.n_sim_steps-1) .* sim.time_step

    power_matrix        = reshape(md.PowerGen, wf.nT, sim.n_sim_steps)
    total_power         = vec(sum(power_matrix, dims=1))
    total_power_trimmed = total_power[(num_timesteps_to_skip+1):end]

    plt_yaw_lim = Plots.plot(
        time_axis, ones(length(time_axis)),
        xlabel         = "Time (s)",
        ylabel         = "Normalised Power\n(-)",
        title          = "Baseline-normalised Farm Power vs Time",
        label          = "Baseline",
        lw             = 2,
        color          = :black,
        legend         = :bottomright,
        legendfontsize = 10,
        tickfontsize = 12,
        guidefontsize = 12,
        grid           = true,
        gridalpha      = 0.4
    )

    yaw_limits = [10.0 15.0 20.0 Inf] # fill in your yaw limit values here

    for yaw_lim in yaw_limits
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
            label = "(Yaw limit = $(yaw_lim)°)",
            lw    = 2,
            alpha = 0.8,
            tickfontsize=12,
            guidefontsize=12,
            legendfontsize=10,
        )
    end

    Plots.plot!(
        plt_yaw_lim,
        size          = (1000, 350),
        left_margin   = 7Plots.mm,
        right_margin  = 15Plots.mm,
        bottom_margin = 5Plots.mm,
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







lambda_curves()
sparsity_curves()
yaw_limit_curves()




