using Timers
using Infiltrator
using ControlPlots
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end
using MAT
using Plots


rerun_sim = false
if rerun_sim




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
#=
yaw_matrix_raw = matread("src/optimisation/result_generation/validation_data/54000_57600h/JHTDB_yaw_data.mat")["yaw_matrix"]
yaw_matrix_raw[:,1] = 3600:7200 # oops 
power_matrix_JHTDB = matread("src/optimisation/result_generation/validation_data/54000_57600h/JHTDB_power_data.mat")["power_matrix"]
=#



yaw_matrix_raw = matread("src/optimisation/result_generation/validation_data/cnbl/CNBL_yaw_data.mat")["yaw_matrix"]
power_matrix_JHTDB = matread("src/optimisation/result_generation/validation_data/cnbl/CNBL_power_data.mat")["power_matrix"] 




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
@time plot_flow_field(wf, X, Y, Z, vis; msr=AddedTurbulence, plt)  #VelReduction, AddedTurbulence, EffWind


end







#=
#SBL CODE










using Statistics

# ── Reshape ──────────────────────────────────────────────────────────────────
n_turbines = 8
n_time     = length(md.PowerGen) ÷ n_turbines
power_matrix_floridyn = [3600:sim.time_step:7200 reshape(md.PowerGen, n_turbines, n_time)']
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




# ── Normalisation ─────────────────────────────────────────────────────────

power_avg_f_normalised = power_avg_f ./ 0.5(power_avg_f[:,1] .+ power_avg_f[:,2])
power_avg_j_normalised = power_avg_j ./ 0.5(power_avg_j[:,1] .+ power_avg_j[:,2])

# ── 8 Subplots ────────────────────────────────────────────────────────────────
window_label = avg_window_sec >= 60 ? "$(avg_window_sec ÷ 60) min" : "$(avg_window_sec) sec"




#power_avg_f_shear_uncorrected = power_avg_f





subplots = map(1:n_turbines) do t
    Plots.plot(
        time_f, power_avg_f[:, t],
        label="FLORIDyn", lw=1.5, ls=:solid, color=4,
        title="Turbine $t",
        ylims=(0, maximum([maximum(power_avg_f), maximum(power_avg_j)])+0.5),  # adjust as needed
        xlabel= "Time (min)",   # only bottom row
        ylabel= mod(t, 2) == 1 ? "Turbine power (MW)" : "",  # only left column
        legend= t == 1 ? :bottomleft : false,
        grid=true, titlefontsize=10,
    )
    #Plots.plot!(time_f, power_avg_f_shear_uncorrected[:, t], label="FLORIDyn (uncorrected shear)", lw=1.5, ls=:solid, color=1)
    Plots.plot!(time_j, power_avg_j[:, t], label="JHTDB LES", lw=1.5, ls=:solid, color=2)
end



p = Plots.plot(subplots...,
               layout=(4, 2),
               size=(900, 1000),
               plot_title="Power Comparison per Turbine, $(window_label) avg",
               leftmargin=5Plots.mm, rightmargin=5Plots.mm, topmargin=1Plots.mm, bottommargin=5Plots.mm)
display(p)









subplots = map(1:n_turbines) do t
    Plots.plot(
        time_f, power_avg_f_normalised[:, t],
        label="FLORIDyn", lw=1.5, ls=:solid,
        title="Turbine $t",
        ylims=(0.0, 1.8),  # adjust as needed
        xlabel= "Time (min)",   # only bottom row
        ylabel= mod(t, 2) == 1 ? "Normalised Power" : "",  # only left column
        legend= t == 1 ? :bottomright : false,
        grid=true, titlefontsize=10,
    )
    Plots.plot!(time_j, power_avg_j_normalised[:, t], label="JHTDB LES", lw=1.5, ls=:solid)
end

p = Plots.plot(subplots...,
               layout=(4, 2),
               size=(900, 1000),
               plot_title="Power Comparison per Turbine ($(window_label) avg)",
               leftmargin=5Plots.mm, rightmargin=5Plots.mm, topmargin=10Plots.mm, bottommargin=5Plots.mm)
display(p)








#Plots.savefig(p, "output/JHTDB_comparison_per_turbine.pdf")




=#





#FOR CNBL

using Statistics

# ── Reshape ──────────────────────────────────────────────────────────────────
n_turbines  = 60
n_rows      = 6
n_cols      = 10
n_timesteps = length(md.PowerGen) ÷ n_turbines
dt          = 1.0   # seconds per timestep

# Reshape power output: [turbine × timestep]
P = reshape(md.PowerGen, n_turbines, n_timesteps)

# ── Column averages ───────────────────────────────────────────────────────────
# Each column k contains row turbines (k-1)*n_rows+1 : k*n_rows
P_col = zeros(n_cols, n_timesteps)
for k in 1:n_cols
    idx = (k-1)*n_rows+1 : k*n_rows
    P_col[k, :] = mean(P[idx, :], dims=1)[:]
end

# ── JHTDB column averages ─────────────────────────────────────────────────────
t_jhtdb_raw = power_matrix_JHTDB[:, 1]
P_jhtdb_raw = power_matrix_JHTDB[:, 2:end] ./ 1e6   # W → MW

P_jhtdb_col = zeros(n_cols, size(P_jhtdb_raw, 1))
for k in 1:n_cols
    idx = (k-1)*n_rows+1 : k*n_rows
    P_jhtdb_col[k, :] = mean(P_jhtdb_raw[:, idx], dims=2)[:]
end

# ── Cut off initialisation period ────────────────────────────────────────────
seconds_to_skip = 600   # 10 min

P_col_trimmed       = P_col[:,      seconds_to_skip+1:end]
P_jhtdb_col_trimmed = P_jhtdb_col[:, seconds_to_skip+1:end]
t_raw_trimmed       = (seconds_to_skip : n_timesteps-1) .* dt          # seconds
t_jhtdb_trimmed     = t_jhtdb_raw[seconds_to_skip+1:end]               # seconds

# ── Averaging window ─────────────────────────────────────────────────────────
avg_window_sec = 5


function window_average(mat, window)
    n          = size(mat, 2)
    n_windows  = n ÷ window
    trimmed    = mat[:, 1:n_windows*window]
    averaged   = [mean(trimmed[k, (w-1)*window+1 : w*window])
                  for k in 1:size(mat,1), w in 1:n_windows]
    return averaged, n_windows
end

P_avg_f, nw_f = window_average(P_col_trimmed,       avg_window_sec)
P_avg_j, nw_j = window_average(P_jhtdb_col_trimmed, avg_window_sec)

# Time axes in minutes (window centres)
time_f = ((1:nw_f) .* avg_window_sec .+ seconds_to_skip) ./ 60
time_j = ((1:nw_j) .* avg_window_sec .+ seconds_to_skip) ./ 60







#choose: absolute or normalised 


# ── 10 subplots ──────────────────────────────────────────────────────────────
window_label = avg_window_sec >= 60 ? "$(avg_window_sec ÷ 60) min" : "$(avg_window_sec) sec"

ylim_max = maximum([maximum(P_avg_f), maximum(P_avg_j)]) + 0.2

subplots = map(1:n_cols) do k
    Plots.plot(
        time_f, P_avg_f[k, :],
        label="FLORIDyn", lw=1.5, ls=:solid, color=1,
        title="Column $k",
        xlims=(10, 60),
        ylims=(0, ylim_max),
        xlabel="Time (min)",
        ylabel=mod(k, 2) == 1 ? "Turbine power (MW)" : "",
        legend=k == 1 ? :bottomright : false,
        grid=true, titlefontsize=10,
    )
    Plots.plot!(time_j, P_avg_j[k, :], label="JHTDB LES", lw=1.5, ls=:solid, color=2)
end

p = Plots.plot(subplots...,
               layout=(5, 2),
               size=(900, 1200),
               plot_title="Column-averaged power comparison, $(window_label) avg",
               leftmargin=5Plots.mm, rightmargin=5Plots.mm,
               topmargin=1Plots.mm, bottommargin=5Plots.mm)
display(p)

Plots.savefig(p, "output/CNBL_comparison_per_column.pdf")


#=



# ── 10 subplots, normalised by Column 1 ──────────────────────────────────────

# Normalise each column by the time-averaged power of Column 1
norm_f = mean(P_avg_f[1, :])
norm_j = mean(P_avg_j[1, :])

P_avg_f_norm = P_avg_f ./ norm_f
P_avg_j_norm = P_avg_j ./ norm_j

window_label = avg_window_sec >= 60 ? "$(avg_window_sec ÷ 60) min" : "$(avg_window_sec) sec"

ylim_max = maximum([maximum(P_avg_f_norm), maximum(P_avg_j_norm)]) + 0.05

subplots = map(1:n_cols) do k
    Plots.plot(
        time_f, P_avg_f_norm[k, :],
        label="FLORIDyn", lw=1.5, ls=:solid, color=1,
        title="Column $k",
        xlims=(10, 60),
        ylims=(0, ylim_max),
        xlabel="Time (min)",
        ylabel=mod(k, 2) == 1 ? "Normalised power (–)" : "",
        legend=k == 1 ? :bottomright : false,
        grid=true, titlefontsize=10,
    )
    Plots.plot!(time_j, P_avg_j_norm[k, :], label="JHTDB LES", lw=1.5, ls=:solid, color=2)
    Plots.hline!([1.0], lw=0.8, ls=:dot, color=:black, label=false)  # reference line at 1
end

p = Plots.plot(subplots...,
               layout=(5, 2),
               size=(900, 1200),
               plot_title="Normalised column power (÷ Col 1 mean), $(window_label) avg",
               leftmargin=5Plots.mm, rightmargin=5Plots.mm,
               topmargin=1Plots.mm, bottommargin=5Plots.mm)
display(p)

Plots.savefig(p, "output/CNBL_comparison_per_column_normalised.pdf")


=#


nothing






