#this file runs the optimisation several times to estimate the variance of the optimiser solution:

using JLD2
using Plots
using Statistics

num_runs = 10;
all_results = []

re_run = false

if re_run

    for i in 1:num_runs
        local result, all_traces, total_time = optimise_yaws();
        push!(all_results, (result, all_traces, total_time));
    end
end

#estimate variance off all result.minimum

min_vals = [r[1].minimum for r in all_results]

println("Mean:   ", mean(min_vals))
println("Var:    ", var(min_vals))
println("Std:    ", std(min_vals))


