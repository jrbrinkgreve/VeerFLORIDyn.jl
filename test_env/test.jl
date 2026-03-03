# the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)


#set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();
set.induction_mode = Induction_TGC();

con.induction = "TGC"


#AAAAAAAAA note: if we keep a between 0.2 and 0.33, we can approximate with constant lambda of 7.6
#meaning we can do axial induction control in this range without changing the code dynamics!!!
# == TO DO FOR LATER

#run initial conditions
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
wf = initSimulation(wf, sim);
con.yaw_data = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate yaw matrix with zeros to prevent error on 'nothing'
con.induction_data = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate induction matrix with zeros to prevent error on 'nothing'


#disable online visualisation
vis.online = false 


#plot flowfield
#construct_yaw_matrix_dynamic!(con.yaw_data, result.minimizer, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
construct_yaw_matrix_dynamic!(con.yaw_data, generate_initial_guess(sim, wind, wf, set_num_yaw_changes), sim, wf, set_num_yaw_changes, set_max_yaw_rate)






construct_axial_induction_matrix_dynamic!(con.induction_data, result.minimizer, sim, wf, set_num_a_changes)
#construct_axial_induction_matrix_dynamic!(con.induction_data, generate_initial_guess_aic(sim, wind, wf, set_num_a_changes), sim, wf, set_num_a_changes)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
#Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
#plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)
include("plot_power.jl")
