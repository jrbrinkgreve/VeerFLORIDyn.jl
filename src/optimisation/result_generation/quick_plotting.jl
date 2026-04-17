#=

#region fold initialisation and setup
#close old plots:

#get the visualisation and settings
vis_file = "data/vis_default.yaml"
settings_file = "data/REALWF_CONTROLTEST_VEER.yaml"   #custom data file with veer specification
vis = Vis(vis_file)
plt=nothing


# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
con.yaw_data = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate yaw matrix with zeros to prevent error on 'nothing'
#run initial conditions
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 



#plot flowfield
construct_yaw_matrix_dynamic!(con.yaw_data, x0, sim, wf, opt_set)
#construct_yaw_matrix_dynamic!(con.yaw_data, result.minimizer, sim, wf, opt_set)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

#top-view
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)






=#
using CSV
using DataFrames
using Dates
using Plots
using Plots.PlotMeasures




const PLOT_FILE  = joinpath(@__DIR__, "../../../data/2021_9T_Data/gains_plot.pdf")  # PDF for publication

function plot_results(wind_directions, baseline_power_avgs, optimized_power_avgs, energy_increases)
    xticks = (wind_directions, string.(wind_directions) .* "°")
 
    p1 = plot(
        wind_directions, [baseline_power_avgs, optimized_power_avgs],
        label       = ["Baseline" "Optimized"],
        xlabel      = "Wind Direction",
        ylabel      = "Average Power (W)",
        marker      = :circle,
        linewidth   = 2,
        xticks      = xticks,
        legend      = :topright,
        dpi         = 300,
    )
 
    p2 = bar(
        wind_directions, energy_increases,
        xlabel    = "Wind Direction",
        ylabel    = "Energy Gain (%)",
        xticks    = xticks,
        legend    = false,
        dpi       = 300,
    )
 
    p = plot(p1, p2, layout = (2, 1), size = (600, 700), dpi = 300)
    savefig(p, PLOT_FILE)
    println("Plot saved to '$PLOT_FILE'")
end
 
plot_results(wind_directions, baseline_power_avgs, optimized_power_avgs, energy_increases)