#this file functions as a recurring call to the yaw controller for a longer MPC formulation

#setup code
hours_to_simulate = 24






#
settings_file = "data/MPC_formulation.yaml"   #custom data file with veer specification



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





















































