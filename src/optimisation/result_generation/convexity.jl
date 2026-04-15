#analyse the convexity of cost_func(x) around x0



using LinearAlgebra

"""
    finite_diff_hessian(f, x0; h=1e-5)

Estimate the Hessian of `f` at `x0` using central finite differences.
Requires n*(n+1)/2 + 1 evaluations (symmetric Hessian), where n = length(x0).

Uses the formula:
  H[i,j] = (f(x+eᵢh+eⱼh) - f(x+eᵢh) - f(x+eⱼh) + f(x)) / h²  for i ≠ j
  H[i,i] = (f(x+eᵢh) - 2f(x) + f(x-eᵢh)) / h²                  for diagonal
"""
function finite_diff_hessian(f, x0; h=1e-5)
    n = length(x0)
    H = zeros(n, n)
    f0 = f(x0)

    # Diagonal terms: second-order central difference
    for i in 1:n
        xp = copy(x0); xp[i] += h
        xm = copy(x0); xm[i] -= h
        H[i, i] = (f(xp) - 2f0 + f(xm)) / h^2
    end

    # Off-diagonal terms: mixed partial via cross differences
    for i in 1:n, j in (i+1):n
        xpp = copy(x0); xpp[i] += h; xpp[j] += h
        xpm = copy(x0); xpm[i] += h; xpm[j] -= h
        xmp = copy(x0); xmp[i] -= h; xmp[j] += h
        xmm = copy(x0); xmm[i] -= h; xmm[j] -= h
        H[i, j] = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4h^2)
        H[j, i] = H[i, j]  # symmetry
    end

    return H
end

"""
    check_convexity(f, x0; h=1e-5, verbose=true)

Estimate the Hessian of `f` at `x0` and check for positive semi-definiteness
as a local convexity indicator.

Returns a NamedTuple with:
  - `eigenvalues`    : sorted eigenvalues of the estimated Hessian
  - `lambda_min`     : smallest eigenvalue
  - `lambda_max`     : largest eigenvalue
  - `is_psd`         : true if all eigenvalues ≥ -tol (locally convex)
  - `hessian`        : the estimated Hessian matrix
  - `n_evaluations`  : number of function evaluations used
"""
function check_convexity(f, x0; h=1e-3, tol=1e-8, verbose=true)
    n = length(x0)
    n_evals = 1 + n + 2 * (n * (n-1) ÷ 2)  # 1 (f0) + n (diag) + 2*(off-diag pairs)

    verbose && println("=" ^ 55)
    verbose && println("  Convexity Check via Hessian Eigenvalue Analysis")
    verbose && println("=" ^ 55)
    verbose && println("  Dimension  : $n")
    verbose && println("  Step size  : h = $h")
    verbose && println("  Evaluations: $n_evals")
    verbose && println("  Est. time  : $(round(n_evals * 0.1, digits=1))s at 0.1s/eval")
    verbose && println("-" ^ 55)

    H = finite_diff_hessian(f, x0; h=h)

    # Symmetrize to counteract floating-point asymmetry
    H_sym = (H + H') / 2

    eigs = eigvals(H_sym)
    sort!(eigs)

    λ_min = eigs[1]
    λ_max = eigs[end]
    is_psd = λ_min ≥ -tol

    if verbose
        println("  λ_min      : $(round(λ_min, sigdigits=6))")
        println("  λ_max      : $(round(λ_max, sigdigits=6))")
        println("  Condition # : $(round(abs(λ_max / λ_min), sigdigits=4))")
        println("-" ^ 55)

        if is_psd
            println("  ✓ Locally CONVEX at x0 (Hessian is PSD)")
        else
            n_neg = count(e -> e < -tol, eigs)
            println("  ✗ NOT locally convex at x0")
            println("    $n_neg negative eigenvalue(s) found")
            println("    Most negative: $(round(λ_min, sigdigits=6))")
        end
        println("=" ^ 55)
    end

    return (
        eigenvalues  = eigs,
        lambda_min   = λ_min,
        lambda_max   = λ_max,
        is_psd       = is_psd,
        hessian      = H_sym,
        n_evaluations = n_evals
    )
end


# ---------------------------------------------------------------
# Example usage — replace cost_func and x0 with your own
# ---------------------------------------------------------------

# A simple 40D convex function for demonstration


sol = check_convexity(cost_func, result.minimizer)

# Access results programmatically
println("\nFull eigenvalue spectrum (sorted):")
println(round.(sol.eigenvalues, sigdigits=4))