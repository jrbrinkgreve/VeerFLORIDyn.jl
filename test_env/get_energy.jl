println("-----------")
println("Optimisation results:")
println(result.minimizer)
println()
println(result.minimum)



# the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

#create settings struct
#
#set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set = Settings(wind, sim, con, false, false)

set.enable_veer = true
set.control_mode = Yaw_Optimisation();


wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);

#run initial conditions
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 

yaw_matrix = construct_yaw_matrix_dynamic(result.minimizer, sim, wf, set_num_yaw_changes)
#plotting:
con.yaw_data = yaw_matrix
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
#Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
#plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)


println(-sum(md.PowerGen))
println()
println("These numbers mismatching is not right!")
println("-----------")
