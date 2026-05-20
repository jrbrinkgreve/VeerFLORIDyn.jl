# code for plotting the cost function:
# three turbines, plot angles of first 2
# slightly mis-positioned with respect to the wind direction
# can also vary veer magnitude for cost function changes

using Infiltrator
using Plots
using JLD2
using Base.Threads


include("../functions.jl")
include("../../../examples/remote_plotting.jl")


# reset states:
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation()

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf = initSimulation(wf, sim)

x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)
con.yaw_data = construct_yaw_matrix_dynamic(x0, sim, wf, set_num_yaw_changes, set_max_yaw_rate)





# Pure mutation, no captured state needed — correct as-is
function construct_yaw_matrix_static!(buffer, x1, x2)
    buffer[:, 2] .= x1
    buffer[:, 3] .= x2
    return nothing
end




# All dependencies are now explicit parameters — no globals captured
function calc_power(md, nT::Int, n_sim_steps::Int, num_timesteps_to_skip)

    return sum(md.PowerGen[(nT*num_timesteps_to_skip+1):end] ) / (nT * (n_sim_steps - num_timesteps_to_skip)) * 1000.0


end

    
function run_cost_sweep(set, wf, wind, sim, con, floridyn, floris, plt, vis, x_array, y_array)

    n_x     = length(x_array)
    n_y     = length(y_array)
    matrix  = zeros(n_x, n_y)
    total   = n_x * n_y

    for (i, x1) in enumerate(x_array)

        for (j, x2) in enumerate(y_array)
            construct_yaw_matrix_static!(con.yaw_data, x1, x2)
            wf_local, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
            #matrix[i, j] = calc_power(md, wf_local.nT, sim.n_sim_steps)

            done = (i - 1) * n_y + j
            @info "[$done/$total] x1=$(round(x1, digits=1))°  x2=$(round(x2, digits=1))°  power=$(round(matrix[i,j], digits=2)) kW"

        end
    end

    return matrix
end


function run_cost_sweep_multithreaded(set, wf, wind, sim, con, floridyn, floris, plt, vis, x_array, y_array)
    n_x    = length(x_array)
    n_y    = length(y_array)
    matrix = zeros(n_x, n_y)
    total  = n_x * n_y

    done_count = Atomic{Int}(0)
    indices    = [(i, j) for i in 1:n_x for j in 1:n_y]

    @threads for idx in indices
        i, j = idx

        # Fresh copy per iteration — no state leaks between tasks
        local_con      = deepcopy(con)
        local_floridyn = deepcopy(floridyn)
        local_wf       = deepcopy(wf)

        x1 = x_array[i]
        x2 = y_array[j]

        construct_yaw_matrix_static!(local_con.yaw_data, x1, x2)

        wf_result, md, mi = run_floridyn(
            plt, set, local_wf, wind, sim, local_con, vis, local_floridyn, floris
        )

        matrix[i, j] = calc_power(md, wf_result.nT, sim.n_sim_steps, 125)

        n_done = atomic_add!(done_count, 1) + 1
        @info "[$n_done/$total] thread=$(threadid())  x1=$(round(x1,digits=1))°  x2=$(round(x2,digits=1))°  power=$(round(matrix[i,j],digits=2)) kW"
    end

    return matrix
end


N = 161
wind_dir_general = 195.0
sweep_size = 80.0

x_array = range(wind_dir_general - sweep_size, wind_dir_general + sweep_size, length=N) #assume 1 wind direction
y_array = range(wind_dir_general - sweep_size, wind_dir_general + sweep_size, length=N)

#comment out to not rerun
#matrix = run_cost_sweep_multithreaded(set, wf, wind, sim, con, floridyn, floris, plt, vis, x_array, y_array)




#=


#plot flowfield
construct_yaw_matrix_dynamic!(con.yaw_data, x0, sim, wf, opt_set)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

#top-view
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)
=#


# wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
# baseline_power_avg = calc_power(md, wf.nT, sim.n_sim_steps)




function plot_cost_landscape(matrix, x_array, y_array, wind)
    x_labels = collect(x_array .- wind.dir[1,2])
    y_labels = collect(y_array .- wind.dir[1,2])

    # --- Contour plot ---
    p1 = contourf(
        x_labels, y_labels, matrix',
        xlabel       = "Turbine 1 Yaw Misalignment (°)",
        ylabel       = "Turbine 2 Yaw Misalignment (°)",
        title        = "Power Landscape",
        color        = :inferno,
        levels       = 25,
        aspect_ratio = :equal,
        grid         = true,
        xlims        = (minimum(x_labels), maximum(x_labels)),
        ylims        = (minimum(y_labels), maximum(y_labels)),
        clims        = (3.0, 4.0),
    )

    # --- Mark the maximum (best operating point) ---
    best_idx = argmax(matrix)
    best_x   = x_labels[best_idx[1]]
    best_y   = y_labels[best_idx[2]]

    scatter!(p1, [best_x], [best_y],
        marker            = :star5,
        markersize        = 10,
        color             = :black,
        markerstrokecolor = :black,
        markerstrokewidth = 1,
        label             = "Max power ($(round(best_x,digits=1))°, $(round(best_y,digits=1))°)",
        legend            = :best,
        legendfontsize    = 8,
    )

    scatter!(p1, [0],[0],
        marker            = :square,
        markersize        = 5,
        color             = :black,
        markerstrokecolor = :black,
        markerstrokewidth = 1,
        label             = "Baseline (0.0°, 0.0°)",
        legend            = :best,
        legendfontsize    = 8,
    )

    # --- Colorbar label (workaround for gr() overlap bug) ---
    annotate!(p1, maximum(x_labels) * 1.40, 0,
        text("Turbines Avg Power (MW)", 9, :center, rotation=90))

    p = plot(p1; size=(600, 500), margin=5Plots.mm)
    display(p)
    savefig(p, "cost_landscape.pdf")
    return p
end

#plot_cost_landscape(matrix*1e-3, x_array, y_array, wind)




using Printf # Required for formatting ticks

function plot_cost_landscape_relative(matrix, x_array, y_array, wind, veer)
    x_labels = collect(x_array .- wind.dir[1,2])
    y_labels = collect(y_array .- wind.dir[1,2])

    # 1. Normalize
    base_idx_x = argmin(abs.(x_labels))
    base_idx_y = argmin(abs.(y_labels))
    baseline_power = matrix[base_idx_x, base_idx_y]
    relative_matrix_pct = ((matrix .- baseline_power) ./ baseline_power) .* 100

    # 2. Get Max point data
    best_idx = argmax(relative_matrix_pct)
    best_x   = x_labels[best_idx[1]]
    best_y   = y_labels[best_idx[2]]
    max_gain = relative_matrix_pct[best_idx]

    # --- Contour plot ---
    p1 = contourf(
        x_labels, y_labels, relative_matrix_pct',
        xlabel       = "Turbine 1 Yaw Misalignment (°)",
        ylabel       = "Turbine 2 Yaw Misalignment (°)",
        title        = "Relative Power Gain at $(veer)°/m Veer",
        color        = :balance,
        levels       = 25,
        aspect_ratio = :equal,
        xlims        = (minimum(x_labels), maximum(x_labels)),
        ylims        = (minimum(y_labels), maximum(y_labels)),
        clims        = (-10, 10), # Your tweaked limits
        colorbar_title = "\nRelative Power Change (%)", # Adds title to colorbar directly
        colorbar_formatter = y -> @sprintf("%.1f%%", y)
    )

    # --- Mark the maximum with coordinates in the label ---
    # Using @sprintf to keep the label clean
    max_label = @sprintf("Max Gain: %.1f%% at (%.1f°, %.1f°)", max_gain, best_x, best_y)

    scatter!(p1, [best_x], [best_y],
        marker            = :star5,
        markersize        = 10,
        color             = :black,
        label             = max_label
    )

    # --- Baseline ---
    scatter!(p1, [0.0], [0.0],
        marker            = :square,
        markersize        = 5,
        color             = :black,
        markerstrokecolor = :black,
        label             = "Baseline (0.0°, 0.0°)"
    )

    # no-veer optimum:
    scatter!(p1, [22.0], [23.0],
        marker            = :diamond,
        markersize        = 8,
        color             = :white,
        markerstrokecolor = :black,
        markerstrokewidth = 1,
        label             = "Optimum at 0°/m Veer (22.0°, 23.0°)"
    )

    # --- Final Layout Adjustment ---
    # Moving legend to :topleft avoids the colorbar on the right entirely
    # Alternatively, use :outertopright but increase the right margin
    p = plot(p1;  
        size=(600, 500), 
        margin=2Plots.mm, 
        legend=:topleft, # Change to :topleft to avoid the right side clash
        legendfontsize=8
    )
    
    display(p)
    savefig(p, "output/$(veer)veer_cost_landscape_relative.pdf")
    return p
end

using JLD2



plot_cost_landscape_relative(matrix, x_array, y_array, wind, floris.veer_gradient)