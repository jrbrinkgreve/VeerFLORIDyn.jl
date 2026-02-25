# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using Timers
using FLORIDyn, TerminalPager, DistributedNext
if Threads.nthreads() == 1; using ControlPlots; end

#_ , vis_file = get_default_project()[2:3]
#vis_file = "data/vis_54T.yaml"
vis_file = "data/vis_54T.yaml"
settings_file = "data/REALWF_54T_CONTROLTEST_VEER.yaml"   #custom data file with veer specification
vis = Vis(vis_file)
vis.show_plots = false  # Enable/disable showing plots during simulation
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end

# Automatic parallel/threading setup
tic()
include("../examples/remote_plotting.jl")
toc()


include("functions.jl")


# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)


# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();


x = result.minimizer




wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
toc()

con.yaw_data = construct_yaw_matrix_dynamic(x, sim, wf, set_num_yaw_changes, set_max_yaw_rate)

vis.online = true
# Clean up any existing PNG files in video folder before starting
cleanup_video_folder()
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
nothing

# was: 61.696510 seconds (7.74 G allocations: 619.266 GiB, 29.02% gc time, 6 lock conflicts, 1.70% compilation time)
# now: 37.828030 seconds (831.03 M allocations: 101.482 GiB, 11.67% gc time, 6 lock conflicts, 2.96% compilation time)