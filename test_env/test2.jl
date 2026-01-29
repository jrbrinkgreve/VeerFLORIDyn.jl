using Evolutionary
using .Threads

# Define objective function
function rosenbrock(x)
    #test parallel calls to floridyn: params already init
    run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris);
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

# Starting point
x0 = [0.0, 0.0]

# Run optimization
# Create custom CMAES instance

opts = Evolutionary.Options(
    iterations = 1000,
    abstol = 1e-8,
    reltol = 1e-8,
    show_trace = true,
    store_trace = true,
    parallelization = :thread
)

result = Evolutionary.optimize(rosenbrock, x0, CMAES(), opts)



# Extract results
Evolutionary.minimizer(result)  # Optimal point
Evolutionary.minimum(result)    # Optimal value
