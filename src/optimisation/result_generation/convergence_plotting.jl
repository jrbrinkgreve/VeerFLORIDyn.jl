function plot_optimiser_trace(; baseline=0.0)
    threshold = 0
    p = Plots.plot(
        title="CMA-ES Convergence",
        xlabel="Iteration",
        ylabel="-P_avg (W)",
        legend=:topright,
        background_color_legend=:white,
        foreground_color_legend=:black,
        tickfontsize=12,
        guidefontsize=12,
        legendfontsize=12,
        left_margin=10Plots.mm,
        bottom_margin=5Plots.mm,
        right_margin=5Plots.mm,
        size=(700, 350)
    )

    # Baseline horizontal line
    Plots.hline!(p, [baseline],
        label="Baseline",
        color=:gray,
        linewidth=2,
        linestyle=:dash)

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
    savefig(p, "output/convergence_plot.pdf")
end

function plot_distance_to_optimum(; baseline=0.0, optimum = 0.0)

    # 90% of the way from baseline down to optimum
    target_down = baseline + 0.9 * (optimum - baseline)
    # Same distance above baseline
    target_up   = baseline - 0.9 * (optimum - baseline)

    p = Plots.plot(
        #title="Distance to Optimum",
        xlabel="Iteration",
        ylabel="-P (kW)",
        legend=:topright,
        background_color_legend=:white,
        foreground_color_legend=:black,
        tickfontsize=12,
        guidefontsize=12,
        legendfontsize=12,
        left_margin=10Plots.mm,
        bottom_margin=5Plots.mm,
        right_margin=5Plots.mm,
        size=(700, 350),
        ylims=(baseline - 1.1*abs(baseline - optimum), baseline + 1.1*abs(baseline - optimum))
    )

    # Baseline and optimum reference lines
    Plots.hline!(p, [baseline],
        label="Baseline",
        color=:gray,
        linewidth=4,
        linestyle=:dash)
    Plots.hline!(p, [optimum],
        label="Global optimum",
        color=:royalblue,
        linewidth=4,
        linestyle=:dash)

    # 90% threshold below baseline (toward optimum)
    Plots.hline!(p, [target_down],
        label="90% to optimum",
        color=:darkorange,
        linewidth=4,
        linestyle=:dot)
                
        let offset = 0
            for (i, trace) in enumerate(all_traces)
                iters  = [t.iteration + offset for t in trace]
                values = [t.value              for t in trace]
                Plots.plot!(p, iters, values, lw=4, color=:black,
                    label = i == 1 ? "CMA-ES with restarts" : false)
                offset = isempty(iters) ? offset : iters[end]
            end
        end

    display(p)
    savefig(p, "output/distance_to_optimum.pdf")
end

# Call with your known baseline, e.g.:
#plot_optimiser_trace(baseline=-3916.40)
plot_distance_to_optimum(baseline=-3916.40, optimum=-3983.64)


nothing