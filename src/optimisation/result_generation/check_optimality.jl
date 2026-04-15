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

# ── core routine ────────────────────────────────────────────────────────────

"""
    slice_1d(f, x★, dim; num_points, half_width)

Return `(ts, vals)` where `ts` are scalar offsets from `x★[dim]` and `vals`
are the corresponding cost-function values along that axis slice.
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
        x[dim]   = x★[dim]          # restore
    end
    return collect(ts), vals
end

"""
    check_optimality(f, x★;
                     num_points = 200,
                     half_width  = 1.0,
                     atol        = 1e-6,
                     plot_slices = true,
                     layout_cols = 3)

For every dimension of `x★`, evaluate `f` along the 1-D slice through `x★`
and check that the minimum of each slice is achieved at (or very near) 0
offset — i.e. at `x★` itself.

Returns a `NamedTuple` with fields:
- `passed`        – `true` if all slices pass the optimality test
- `dim_results`   – per-dimension `NamedTuple`s (slice_min_offset, slice_min_val, optimal_val, passed)
- `fig`           – the Plots figure (or `nothing` if `plot_slices=false`)
"""
function check_optimality(f, x★::AbstractVector;
                          num_points::Int  = 200,
                          half_width::Real = 1.0,
                          atol::Real       = 1e-6,
                          plot_slices::Bool = true,
                          layout_cols::Int  = 6)

    n          = length(x★)
    f★         = f(x★)
    dim_results = Vector{NamedTuple}(undef, n)

    println("="^60)
    println("  1-D Slice Optimality Check")
    println("  Candidate minimizer f(x★) = $f★")
    println("="^60)

    plots_list = plot_slices ? Vector{Plots.Plot}(undef, n) : nothing

    all_passed = true
    for d in 1:n
        ts, vals = slice_1d(f, x★, d; num_points, half_width)

        min_val    = minimum(vals)
        min_offset = ts[argmin(vals)]
        passed_d   = abs(min_offset) ≤ atol || min_val ≥ f★ - atol

        all_passed = all_passed && passed_d
        status     = passed_d ? "✓ PASS" : "✗ FAIL"

        @printf("  dim %3d | f(x★)= %+.6e | slice_min= %+.6e at offset %+.3e | %s\n",
                d, f★, min_val, min_offset, status)

        dim_results[d] = (dim            = d,
                          slice_min_offset = min_offset,
                          slice_min_val    = min_val,
                          optimal_val      = f★,
                          passed           = passed_d)

        if plot_slices
            col = passed_d ? :steelblue : :crimson
            p   = plot(ts, vals;
                       label     = "f along dim $d",
                       color     = col,
                       lw        = 2,
                       xlabel    = "offset from x★[$d]",
                       ylabel    = "cost",
                       title     = "dim $d  ($status)",
                       titlefont = font(9),
                       legend    = false)
            vline!(p, [0.0]; color=:black, lw=1.5, ls=:dash, label="x★")
            hline!(p, [f★];  color=:grey,  lw=1,   ls=:dot,  label="f(x★)")
            plots_list[d] = p
        end
    end

    println("="^60)
    println("  Overall: ", all_passed ? "ALL DIMENSIONS PASS ✓" : "SOME DIMENSIONS FAILED ✗")
    println("="^60)

    fig = nothing
    if plot_slices
        cols = min(layout_cols, n)
        rows = cld(n, cols)
        subplot_w = 250
        subplot_h = 200
        fig  = plot(plots_list...;
                    layout    = (rows, cols),
                    size      = (subplot_w * cols, subplot_h * rows),
                    margin    = 3Plots.mm,
                    titlefont = font(7),
                    guidefont = font(6),
                    tickfont  = font(5))
        display(fig)
    end

    return (passed = all_passed, dim_results = dim_results, fig = fig)
end


# ── example / demo ──────────────────────────────────────────────────────────

# Run the check
report = check_optimality(cost_func, result.minimizer;
                           num_points  = 30,
                           half_width  = 0.08,
                           atol        = 1e-4,
                           plot_slices = true,
                           layout_cols = 8)

# Inspect per-dimension results programmatically
for r in report.dim_results
    println("dim $(r.dim): passed=$(r.passed), slice_min_offset=$(round(r.slice_min_offset, sigdigits=4))")
end



