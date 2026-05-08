#this file runs the optimisation several times to estimate the variance of the optimiser solution:

using JLD2
using Plots
using Statistics

num_runs = 50;

re_run = false
if !@isdefined(all_results)
    all_results = []
end 

if re_run

    for i in 1:num_runs
        println("Run $i / $num_runs")
        include("../parallel_yaw_optimisation.jl")
        push!(all_results, (result, all_traces, total_time));
    end
end

#estimate variance off all result.minimum

min_vals = [r[1].minimum for r in all_results]
t_vals = [r[1].minimizer[1:3] for r in all_results] * 4100
yaw_vals = [r[1].minimizer[4:end] for r in all_results] * 360


p_increase = (-3916.4 .- min_vals) / 3916.4

println("Mean:   ", mean(min_vals))
println("Std:    ", std(min_vals))



println("Mean:   ", mean(t_vals))
println("Std:    ", std(t_vals))



println("Mean:   ", mean(yaw_vals))
println("Std:    ", std(yaw_vals))



println("Mean:   ", mean(p_increase) * 100, " %")
println("Std:    ", std(p_increase) * 100, " %" )


control_instructions = [r[1].minimizer for r in all_results]


nothing

#@save "../optim_data/solution_variance_tests.jld2" all_results min_vals t_vals yaw_vals p_increase
#@load "../optim_data/baseline_init_solution_variance_tests.jld2" all_results min_vals t_vals yaw_vals p_increase
