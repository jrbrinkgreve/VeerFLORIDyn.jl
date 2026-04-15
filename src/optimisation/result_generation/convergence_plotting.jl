
#convergence plot
threshold = 0
p = Plots.plot(title="CMA-ES Convergence", xlabel="Iteration", ylabel="Cost")

let offset = 0
    for (i, trace) in enumerate(all_traces)
        filtered = [(t.iteration, t.value) for t in trace if t.value < threshold]
        iters    = [x[1] + offset for x in filtered]
        values   = [x[2]          for x in filtered]
        Plots.plot!(p, iters, values, lw=2, label="Run $i")
        offset = isempty(iters) ? offset : iters[end]
    end
end
display(p)
savefig(p, "convergence_plot.pdf")
