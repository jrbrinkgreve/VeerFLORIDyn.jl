using Timers
using Evolutionary
using OhMyThreads: TaskLocalValue
using Infiltrator
using Profile
using BenchmarkTools
using LinearAlgebra
using FLORIDyn, TerminalPager, DistributedNext
using Plots
using Printf
if Threads.nthreads() == 1; using ControlPlots; end


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

# (Re)initialise whenever data changes
init_tlv!()





#endregion

#------------------------------------------------------------------------------------------------------------------
#OPTIMISATION PART

#sim 
set_num_yaw_changes = 4         #N
set_max_yaw_misalignment = 25.0 #deg, for penalising large yaw angles in the cost function, for stability and convergence reasons
set_lambda_l1 = 0.0            #1e3  #units: cost PER DEGREE, per turbine, PER SECOND #relative to the beneficial term average kW per turbine #typical value 1e3, can play around with this
set_lambda_l1_hard_limit = Inf  #a limit on the maximum total yaw change in a simulation, in degrees
set_max_yaw_rate = 1.0          #deg/s
set_objective = totalEnergyObjective   #totalEnergyObjective or powerTrackingObjective al
set_num_timesteps_to_skip = 125     #skip the first N timesteps for wake effects to propagate, approx time between wake interactions

#optimiser convergence/ fidelity  
set_cmaes_lambda_multiplier = 4    # 4       #multiplier for the default lambda, which is 4 + 3 * log(N), N is dim of problem
set_num_optimiser_runs = 4        # 4       #number of automatic restarts for CMA-ES
set_iterations = 50                # 50     #number of iterations for CMAES
set_sigma0 = 0.05                  # 0.05   #for time optim, 0.01 works well in second run for yaws!!     # set to 30% of the search range, and for yaw convergence: first 0.05 then 0.01 
set_sigma0_secondary = 0.02        # 0.02   #for second run with yaws, to reduce the search area and converge faster, can play around with this #0.01
set_sigma0_final = 0.01            # 0.01  #for final run 



#initial state
x0 = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)


results = []

for veer in 0.0:0.01:0.1




        
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

    floris.veer_gradient = veer
    
    construct_yaw_matrix_dynamic!(con.yaw_data, x0, sim, wf, opt_set)
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    
    optimized_power_avg = sum(md.PowerGen[(wf.nT*num_timesteps_to_skip+1):end]) / 
                          (wf.nT * (sim.n_sim_steps - num_timesteps_to_skip)) * 1000.0
    
    push!(results, (veer=veer, power=optimized_power_avg))
    @printf("veer = %.4f | P_opt = %.2f kW\n", veer, optimized_power_avg)
end

# Summary
println("\n--- Summary ---")
println("veer\t\tP_opt (kW)")
for r in results
    @printf("%.4f\t\t%.2f\n", r.veer, r.power)
end

