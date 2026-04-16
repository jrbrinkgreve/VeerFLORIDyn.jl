"""
    random_direction_optimality_check.jl

Verify optimality of a candidate minimizer by walking 1D slices of the cost
function along N random directions in the full high-D space.

Instead of axis-aligned slices (O(n) checks), this samples random unit vectors
and checks that x★ is a local minimum along each random cross-section.
"""

using Printf
using Base.Threads
using LinearAlgebra

"""
    slice_random_direction(f, x★, direction; num_points, half_width)

Evaluate f along the line x★ + t * direction for t ∈ [-half_width, half_width].
`direction` should be a unit vector.
"""
function slice_random_direction(f, x★::AbstractVector, direction::AbstractVector;
                                 num_points::Int  = 200,
                                 half_width::Real = 1.0)
    ts   = range(-half_width, half_width; length = num_points)
    vals = Vector{Float64}(undef, num_points)
    for (i, t) in enumerate(ts)
        vals[i] = f(x★ .+ t .* direction)
    end
    return collect(ts), vals
end

"""
    check_convexity_random(f, x★; N, num_points, half_width, atol, seed)

Check that x★ is a local minimum along N random directions in the ambient space.

For each direction d (a random unit vector), we verify that the 1D function
    g(t) = f(x★ + t·d)
has its minimum at t ≈ 0, i.e. x★ is optimal along that slice.

Arguments
---------
- `f`           : cost function R^n → R
- `x★`          : candidate minimizer
- `N`           : number of random directions to sample
- `num_points`  : resolution of each 1D scan
- `half_width`  : scan range is [-half_width, +half_width] around x★
- `atol`        : tolerance for declaring the minimum is at t=0
- `seed`        : RNG seed for reproducibility (nothing = random)
- `n_threads`   : whether to use multithreading

Returns a NamedTuple with fields:
- `passed`       : Bool, true iff all directions passed
- `dir_results`  : Vector of per-direction result NamedTuples
- `directions`   : the sampled unit vectors (n_dims × N matrix)
- `pass_rate`    : fraction of directions that passed
"""
function check_convexity_random(f, x★::AbstractVector;
                                 N::Int           = 100,
                                 num_points::Int  = 200,
                                 half_width::Real = 1.0,
                                 atol::Real       = 1e-6,
                                 seed             = 42,
                                 n_threads::Bool  = true)

    n  = length(x★)
    f★ = f(x★)

    # --- Sample N random unit vectors ---
    rng        = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)
    raw_dirs   = randn(rng, n, N)           # n × N matrix of Gaussian vectors
    directions = mapslices(v -> v ./ norm(v), raw_dirs; dims=1)  # normalise each column

    dir_results = Vector{NamedTuple}(undef, N)

    println("="^72)
    println("  Random-Direction Convexity / Optimality Check")
    @printf("  Dimensions: %d  |  Directions: %d  |  f(x★) = %.6e\n", n, N, f★)
    println("-"^72)

    done_count = Atomic{Int}(0)

    run_check = function(k)
        dir = directions[:, k]

        ts, vals = slice_random_direction(f, x★, dir; num_points, half_width)

        min_val    = minimum(vals)
        min_idx    = argmin(vals)
        min_offset = ts[min_idx]

        # Pass if: the slice minimum is at (or very near) t=0, OR the slice
        # minimum is not meaningfully lower than f★ (flat / numerical noise).
        passed_k = abs(min_offset) ≤ atol || min_val ≥ f★ - atol
        status   = passed_k ? "✓" : "✗ FAIL"

        dir_results[k] = (
            direction_index  = k,
            min_offset       = min_offset,
            min_val          = min_val,
            optimal_val      = f★,
            passed           = passed_k,
            direction        = dir,
        )

        n_done = atomic_add!(done_count, 1) + 1
        if !passed_k || mod(n_done, max(1, N ÷ 10)) == 0  # print every 10% + all failures
            @printf("  [%4d/%4d] slice_min= %+.4e at t= %+.3e | %s\n",
                    n_done, N, min_val, min_offset, status)
        end
    end

    if n_threads
        @threads for k in 1:N
            run_check(k)
        end
    else
        for k in 1:N
            run_check(k)
        end
    end

    # --- Aggregate ---
    all_passed = all(r -> r.passed, dir_results)
    n_passed   = count(r -> r.passed, dir_results)
    pass_rate  = n_passed / N

    failures   = filter(r -> !r.passed, dir_results)

    println("-"^72)
    @printf("  Pass rate: %d / %d  (%.1f%%)\n", n_passed, N, 100 * pass_rate)

    if !isempty(failures)
        println("\n  Failed directions (worst offenders):")
        sorted_failures = sort(failures; by = r -> r.min_val)  # worst (lowest) first
        for r in Iterators.take(sorted_failures, min(10, length(failures)))
            @printf("    dir #%d : slice_min= %+.4e at t= %+.3e  (Δ = %.4e below f★)\n",
                    r.direction_index, r.min_val, r.min_offset, f★ - r.min_val)
        end
    end

    println("  Overall: ", all_passed ? "ALL DIRECTIONS PASS ✓" : "SOME DIRECTIONS FAILED ✗")
    println("="^72)

    return (
        passed      = all_passed,
        pass_rate   = pass_rate,
        dir_results = dir_results,
        directions  = directions,
    )
end


# --- Usage ---
import Random   # needed for MersenneTwister

report = check_convexity_random(cost_func, result.minimizer;
             N           = 100,     # number of random directions
             num_points  = 51,     # points per 1D scan
             half_width  = 0.5,     # scan ± 0.5 around x★ along each direction
             atol        = 1.0,     # tolerance (match your original)
             seed        = 42,      # reproducible direction sampling
             n_threads   = true)

# Inspect failures in detail
failures = filter(r -> !r.passed, report.dir_results)
for r in failures
    println("dir #$(r.direction_index): offset=$(round(r.min_offset, sigdigits=4)), " *
            "Δ=$(round(report.dir_results[1].optimal_val - r.min_val, sigdigits=4))")
end