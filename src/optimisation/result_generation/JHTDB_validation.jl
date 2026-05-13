using Timers
using Infiltrator
using ControlPlots
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end
using MAT
using Plots




#_ , vis_file = get_default_project()[2:3]
vis_file = "data/vis_default.yaml"
settings_file = "data/JHTDB_comparison_turbines.yaml"   #custom data file with veer specification

# Load vis settings from YAML file
vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end

# Automatic parallel/threading setup

include("../../../examples/remote_plotting.jl")


# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)
# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();

#load data
yaw_matrix_raw = matread("src/optimisation/result_generation/validation_data/JHTDB_yaw_data.mat")["yaw_matrix"]
TI_matrix = matread("src/optimisation/result_generation/validation_data/JHTDB_TI_data.mat")["TI_array"]
winddir_matrix = matread("src/optimisation/result_generation/validation_data/JHTDB_winddir_data.mat")["winddir_array"]
windspeed_matrix = matread("src/optimisation/result_generation/validation_data/JHTDB_windspeed_data.mat")["windspeed_array"]




#convert to deg and coordinates
yaw_matrix = zeros(size(yaw_matrix_raw))
yaw_matrix[:,1] = yaw_matrix_raw[:, 1]  
yaw_matrix[:,2:end] = 270.0 .- rad2deg.(yaw_matrix_raw[:, 2:end])
con.yaw_data = yaw_matrix
con.yaw = "Optimisation";





wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)





# Run initial conditions
wf = initSimulation(wf, sim)





# toggle for video file creation
make_video = false

if make_video
    cleanup_video_folder()  # Ensure video folder is clean before starting
end
vis.show_plots = !make_video 
vis.online = make_video #false



@time wf, md, mi = run_floridyn(ControlPlots.plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
@time plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)  #VelReduction, AddedTurbulence, EffWind










using Statistics

# ── Reshape ──────────────────────────────────────────────────────────────────
n_turbines = 8
n_time     = length(md.PowerGen) ÷ n_turbines
power_matrix_floridyn = [3600:sim.time_step:7200 reshape(md.PowerGen, n_turbines, n_time)']
power_matrix_JHTDB = matread("src/optimisation/result_generation/validation_data/JHTDB_power_data.mat")["power_matrix"]
power_matrix_JHTDB[:, 2:end] = power_matrix_JHTDB[:, 2:end] .* 1e-6

# ── Cut off initialisation period ────────────────────────────────────────────
seconds_to_skip = 600
power_matrix_floridyn = power_matrix_floridyn[seconds_to_skip:end, :]
power_matrix_JHTDB    = power_matrix_JHTDB[seconds_to_skip:end, :]
n_time = size(power_matrix_floridyn, 1)

# ── Averaging window ─────────────────────────────────────────────────────────
avg_window_sec = 30

function window_average(mat, window, n_turb)
    n = size(mat, 1)
    n_windows = n ÷ window
    trimmed = mat[1:(n_windows * window), :]
    averaged = [mean(trimmed[(w-1)*window+1 : w*window, t])
                for w in 1:n_windows, t in 1:n_turb]
    return averaged, n_windows
end

power_avg_f, n_windows_f = window_average(power_matrix_floridyn[:, 2:end], avg_window_sec, n_turbines)
power_avg_j, n_windows_j = window_average(power_matrix_JHTDB[:, 2:end],    avg_window_sec, n_turbines)

time_f = (1:n_windows_f) .* avg_window_sec ./ 60 .+ 10
time_j = (1:n_windows_j) .* avg_window_sec ./ 60 .+ 10

# ── 8 Subplots ────────────────────────────────────────────────────────────────
window_label = avg_window_sec >= 60 ? "$(avg_window_sec ÷ 60) min" : "$(avg_window_sec) sec"

subplots = map(1:n_turbines) do t
    Plots.plot(
        time_f, power_avg_f[:, t],
        label="FLORIDyn", lw=1.5, ls=:solid,
        title="Turbine $t",
        ylims=(2.5, 12.5),  # adjust as needed
        xlabel= "Time (min)",   # only bottom row
        ylabel= mod(t, 2) == 1 ? "Power (MW)" : "",  # only left column
        legend= t == 1 ? :bottomright : false,
        grid=true, titlefontsize=10,
    )
    Plots.plot!(time_j, power_avg_j[:, t], label="JHTDB LES", lw=1.5, ls=:solid)
end

p = Plots.plot(subplots...,
               layout=(4, 2),
               size=(900, 1000),
               plot_title="Power Comparison per Turbine ($(window_label) avg)",
               leftmargin=5Plots.mm, rightmargin=5Plots.mm, topmargin=10Plots.mm, bottommargin=5Plots.mm)
display(p)



Plots.savefig(p, "output/JHTDB_comparison_per_turbine.pdf")















nothing





