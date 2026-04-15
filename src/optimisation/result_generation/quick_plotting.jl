

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






