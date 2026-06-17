# code for plotting the cost function:
# three turbines, plot angles of first 2
# slightly mis-positioned with respect to the wind direction
# can also vary veer magnitude for cost function changes

using Infiltrator
using Plots
using JLD2
using Base.Threads
using Printf


include("../functions.jl")
include("../../../examples/remote_plotting.jl")


# reset states:
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation()

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf = initSimulation(wf, sim)



#switching time

#load optimal switching time: 


path   = "src/optimisation/optim_data/minimisers/sparsity_4.jld2"
result = load(path)["result"]
x      = result.minimizer

#offset x[1] and x[2] by 0.25 in 100 steps each to plot cost function landscape 

con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)









# All dependencies are now explicit parameters — no globals captured
function calc_power(md, nT::Int, n_sim_steps::Int, num_timesteps_to_skip)

    return sum(md.PowerGen[(nT*num_timesteps_to_skip+1):end] ) / (nT * (n_sim_steps - num_timesteps_to_skip)) * 1000.0


end



# Pure mutation helper — offsets only the first two elements of x
function construct_x_offset!(x_local, x_base, δ1, δ2)
    x_local .= x_base
    x_local[2] += δ1   # 2 is second to last time
    x_local[3] += δ2  # 24 is front turbine yaw at second to last time
    return nothing
end


function run_switching_cost_sweep_multithreaded(
    set, wf, wind, sim, con, floridyn, floris, plt, vis,
    x_base,                    # the minimiser loaded from JLD2
    d1_array, d2_array,        # offset grids for x[1] and x[2]
    num_timesteps_to_skip = 125
)
    n1     = length(d1_array)
    n2     = length(d2_array)
    matrix = zeros(n1, n2)
    total  = n1 * n2

    done_count = Atomic{Int}(0)
    indices    = [(i, j) for i in 1:n1 for j in 1:n2]

    @threads for idx in indices
        i, j = idx

        # Thread-local deep copies — no shared mutable state
        local_con      = deepcopy(con)
        local_floridyn = deepcopy(floridyn)
        local_wf       = deepcopy(wf)
        local_x        = copy(x_base)

        δ1 = d1_array[i]
        δ2 = d2_array[j]

        construct_x_offset!(local_x, x_base, δ1, δ2)

        # Rebuild the yaw matrix from the perturbed x
        local_con.yaw_data = construct_yaw_matrix_dynamic(
            local_x, sim, local_wf, set_num_yaw_changes, set_max_yaw_rate
        )

            
        #wf_result, md, mi = run_floridyn(
        #    plt, set, local_wf, wind, sim, local_con, vis, local_floridyn, floris
        #)

        
        #matrix[i, j] = calc_power(md, wf_result.nT, sim.n_sim_steps, num_timesteps_to_skip)


        matrix[i, j] = cost_func(local_x)
        n_done = atomic_add!(done_count, 1) + 1
        @info "[$n_done/$total] thread=$(threadid())  δ1=$(round(δ1,digits=3))  δ2=$(round(δ2,digits=3))  power=$(round(matrix[i,j],digits=2)) kW"
    end

    return matrix
end


# --- Setup and run ---
N           = 51
offset_size_1 = 0.3
offset_size_2 =  0.3

d1_array = range(-offset_size_1, offset_size_1, length=N)
d2_array = range(-offset_size_2, offset_size_2, length=N)

switching_matrix = run_switching_cost_sweep_multithreaded(
    set, wf, wind, sim, con, floridyn, floris, plt, vis,
    x, d1_array, d2_array
)







# --- Plotting (maps an N x N matrix directly to ±1025 axes) ---
function plot_switching_landscape(matrix)
    # 1. Dynamically create axes that span exactly from -1025 to 1025
    N = size(matrix, 1)
    d1 = range(-1230, 1230, length=N)
    d2 = range(-1230, 1230, length=N)

    # 2. Find the index closest to 0 for the baseline calculation
    base_idx_1 = argmin(abs.(d1))
    base_idx_2 = argmin(abs.(d2))
    baseline   = matrix[base_idx_1, base_idx_2]
    
    # Calculate relative changes (%)
    rel_tmp = ((matrix .- baseline) ./ baseline) .* 100
    rel = clamp.(rel_tmp, -5.0, 5.0)

    best_idx = argmax(rel)
    best_d1  = d1[best_idx[1]]
    best_d2  = d2[best_idx[2]]
    max_gain = rel[best_idx]

    max_label = @sprintf("Max Gain: %.2f%% at (%.3f, %.3f)", max_gain, best_d1, best_d2)

    # 3. Generate the plot
    p1 = contourf(
        d1, d2, rel',
        xlabel             = "Offset on switching time 3 \n (s)",
        ylabel             = "Offset on switching time 4 \n (s)",
        title              = "Joint Switching Time Cost Landscape",
        color              = :balance,
        levels             = 50,
        #aspect_ratio       = :equal,
        xlims              = (-1230, 1230),
        ylims              = (-1230, 1230),
        clims              = (-5, 5),
    
        colorbar_title     = "\nRelative Power Change (%)",
        colorbar_formatter = y -> @sprintf("%.1f%%", y),
    )

      scatter!(p1, [0.0], [0.0],
        marker            = :star,
        markersize        = 7,
        color             = :black,
        markerstrokecolor = :black,
        label             = "optimum",
    
    )

    p = plot!(p1;
        size         = (600, 500),
        margin       = 5Plots.mm,

        legendfontsize = 12,
        titlefontsize  = 12,
        tickfontsize   = 12,
        guidefontsize   = 12,
        colorbar_titlefontsize = 12,
    )

    display(p)
    savefig(p, "output/switching_cost_landscape.pdf")
    return p
end

# To run it, you just pass your matrix:
plot_switching_landscape(switching_matrix)