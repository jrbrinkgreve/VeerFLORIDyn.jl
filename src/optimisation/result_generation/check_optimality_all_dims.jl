"""
    slice_optimality_check.jl

Verify optimality of a candidate minimizer by walking 1D slices of the cost
function along each coordinate axis.  Complexity is O(n · num_points) instead
of the exponential O(num_points^n) full grid.

Usage
-----
1. Paste / import your `cost_func` and `result` before running, or adjust the
   example at the bottom of this file.
2. Call `check_optimality(cost_func, result.minimizer)`.
"""

using Printf
using Plots

"""
    slice_1d(f, x★, dim; num_points, half_width)
"""
function slice_1d(f, x★::AbstractVector, dim::Int;
                  num_points::Int = 200,
                  half_width::Real = 1.0)
    ts   = range(-half_width, half_width; length = num_points)
    vals = Vector{Float64}(undef, num_points)
    x    = copy(x★)
    for (i, t) in enumerate(ts)
        x[dim]   = x★[dim] + t
        vals[i]  = f(x)
        x[dim]   = x★[dim]          
    end
    return collect(ts), vals
end

"""
    check_optimality(f, x★; ...)
"""
function check_optimality(f, x★::AbstractVector;
                          num_points::Int  = 200,
                          half_width::Real = 1.0,
                          atol::Real       = 1e-6,
                          plot_slices::Bool = true,
                          layout_cols::Int  = 5) # 5 cols is ideal for 10 dims

    n           = length(x★)
    f★          = f(x★)
    dim_results = Vector{NamedTuple}(undef, n)

    # --- Header Tracker ---
    println("="^70)
    println("  1-D Slice Optimality Check")
    println("  Target f(x*) = $f★")
    println("-"^70)

    plots_list = plot_slices ? Vector{Plots.Plot}(undef, n) : nothing

    all_passed = true
    for d in 1:n
        ts, vals = slice_1d(f, x★, d; num_points, half_width)

        min_val    = minimum(vals)
        min_idx    = argmin(vals)
        min_offset = ts[min_idx]
        
        # Logic check: is the lowest point in the slice at the center (0.0)?
        passed_d   = abs(min_offset) ≤ atol || min_val ≥ f★ - atol
        all_passed = all_passed && passed_d
        status     = passed_d ? "✓ PASS" : "✗ FAIL"

        # --- Console Progress Tracker ---
        @printf("  dim %3d | f(x★)= %+.4e | slice_min= %+.4e at offset %+.3e | %s\n",
                d, f★, min_val, min_offset, status)

        dim_results[d] = (dim = d, slice_min_offset = min_offset, 
                          slice_min_val = min_val, optimal_val = f★, passed = passed_d)

        if plot_slices
            line_col = passed_d ? :steelblue : :crimson
            bg_col   = passed_d ? :white : RGBA(1.0, 0.92, 0.92, 1.0) # Light red tint for fails
            
            p = plot(ts, vals;
                       label     = "f(x)",
                       color     = line_col,
                       lw        = 2,
                       title     = "Dim $d ($status)",
                       titlefontsize = 9,
                       background_color_subplot = bg_col,
                       xlabel    = "offset",
                       ylabel    = "cost",
                       guidefontsize = 7,
                       tickfontsize  = 6,
                       legend    = false)

            # Visual markers
            vline!(p, [0.0], color=:orange, ls=:dash, lw=1.2, label="x★")
            hline!(p, [f★], color=:black, ls=:dot, alpha=0.5)
            
            # Highlight the actual local minimum found in the slice
            scatter!(p, [min_offset], [min_val], marker=:circle, ms=2, mc=:white, msc=line_col)

            plots_list[d] = p
        end
    end

    # --- Footer Tracker ---
    println("-"^70)
    println("  Overall: ", all_passed ? "ALL DIMENSIONS PASS ✓" : "SOME DIMENSIONS FAILED ✗")
    println("="^70)

    fig = nothing
    if plot_slices
        cols = min(layout_cols, n)
        rows = cld(n, cols)
        
        fig = plot(plots_list...;
                    layout    = (rows, cols),
                    size      = (250 * cols, 220 * rows),
                    margin    = 4Plots.mm,
                    plot_title = "Optimality Slices",
                    plot_titlefontsize = 12)
        display(fig)
    end

    return (passed = all_passed, dim_results = dim_results, fig = fig)
end


using Base.Threads
function check_optimality_multithreaded(f, x★::AbstractVector;
    num_points::Int  = 200,
    half_width::Real = 1.0,
    atol::Real       = 1e-6,
    plot_slices::Bool = true,
    layout_cols::Int  = 5)

    n           = length(x★)
    f★          = f(x★)
    dim_results = Vector{NamedTuple}(undef, n)

    println("="^70)
    println("  1-D Slice Optimality Check (Multithreaded)")
    println("  Target f(x*) = $f★")
    println("-"^70)

    plots_list = plot_slices ? Vector{Plots.Plot}(undef, n) : nothing

    done_count = Atomic{Int}(0)

    @threads for d in 1:n

        # --- LOCAL COPIES (thread safety!) ---
        x_local = copy(x★)   # deepcopy not needed unless x contains objects

        ts   = range(-half_width, half_width; length = num_points)
        vals = Vector{Float64}(undef, num_points)

        for (i, t) in enumerate(ts)
            x_local[d] = x★[d] + t
            vals[i]    = f(x_local)
            x_local[d] = x★[d]
        end

        min_val    = minimum(vals)
        min_idx    = argmin(vals)
        min_offset = ts[min_idx]

        passed_d = abs(min_offset) ≤ atol || min_val ≥ f★ - atol
        status   = passed_d ? "✓ PASS" : "✗ FAIL"

        dim_results[d] = (
            dim = d,
            slice_min_offset = min_offset,
            slice_min_val = min_val,
            optimal_val = f★,
            passed = passed_d
        )

        # --- Progress tracking ---
        n_done = atomic_add!(done_count, 1) + 1
        @printf("  [%d/%d] dim %3d | f(x★)= %+.4e | slice_min= %+.4e at offset %+.3e | %s \n",
                n_done, n, d, f★, min_val, min_offset, status)

        # --- Plotting (thread-safe: each thread writes its own slot) ---
            if plot_slices
            line_col  = passed_d ? :steelblue : :crimson
            bg_col    = passed_d ? :white : RGBA(1.0, 0.92, 0.92, 1.0)

            # Clip indicator spikes for display only
            vals_plot = replace(v -> v > 1e5 ? NaN : v, vals)
            dot_y     = min_val > 1e5 ? NaN : min_val

            p = plot(ts, vals_plot;
                    label     = "f(x)",
                    color     = line_col,
                    lw        = 2,
                    title     = "Dim $d ($status)",
                    titlefontsize = 9,
                    background_color_subplot = bg_col,
                    xlabel    = "offset",
                    ylabel    = "cost",
                    guidefontsize = 7,
                    tickfontsize  = 6,
                    legend    = false)

            # --- NEW: Highlight Infeasible Regions (vspan) ---
            # We find contiguous blocks of NaNs to draw the red spans
            i = 1
            while i <= length(ts)
                if isnan(vals_plot[i])
                    start_val = ts[i]
                    # Walk until we find a non-NaN value or the end of the array
                    while i <= length(ts) && isnan(vals_plot[i])
                        i += 1
                    end
                    end_val = ts[min(i, length(ts))]
                    
                    # Apply the vertical red tint
                    vspan!(p, [start_val, end_val], 
                        color=:red, alpha=0.15, linecolor=:transparent)
                    
                    # Optional: Add a thick red line specifically on the x-axis
                    y_low, _ = ylims(p)
                    plot!(p, [start_val, end_val], [y_low, y_low], 
                        color=:red, lw=3, primary=false)
                else
                    i += 1
                end
    end
    # -----------------------------------------------

    vline!(p, [0.0], color=:orange, ls=:dash, lw=1.2, label="x★")
    hline!(p, [f★],  color=:black,  ls=:dot,  alpha=0.5)
    scatter!(p, [min_offset], [dot_y], marker=:circle, ms=2, mc=:white, msc=line_col)

    plots_list[d] = p
end
    end

    # --- Final aggregation ---
    all_passed = all(r -> r.passed, dim_results)

    println("-"^70)
    println("  Overall: ", all_passed ? "ALL DIMENSIONS PASS ✓" : "SOME DIMENSIONS FAILED ✗")
    println("="^70)

    fig = nothing
    if plot_slices
        cols = min(layout_cols, n)
        rows = cld(n, cols)

    fig = plot(plots_list...;
        layout    = (rows, cols),
        size      = (250 * cols, 260 * rows),   # slightly taller rows
        margin    = 4Plots.mm,
        top_margin    = 8Plots.mm,              # extra vertical breathing room
        bottom_margin = 8Plots.mm,
        plot_title = "Optimality Slices",
        plot_titlefontsize = 16)

        display(fig)
    end
    savefig(fig, "optimality_slices.pdf")
    return (passed = all_passed, dim_results = dim_results, fig = fig)
end




# Usage
report = check_optimality_multithreaded(cost_func, result.minimizer;
                           num_points  = 51,
                           half_width  = 0.25,   #180deg, 
                           atol        = 1e-2,  #tolerance for optimality, this is 1e-2 kW now
                           plot_slices = true,
                           layout_cols = 6)


for r in report.dim_results
    println("dim $(r.dim): passed=$(r.passed), slice_min_offset=$(round(r.slice_min_offset, sigdigits=4))")
end



