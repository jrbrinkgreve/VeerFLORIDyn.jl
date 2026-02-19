using Timers
using Evolutionary
using OhMyThreads: TaskLocalValue
using Infiltrator
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end


#get the visualisation and settings
_ , vis_file = get_default_project()[2:3]
settings_file = "data/CONTROLTEST_VEER.yaml"   #custom data file with veer specification

vis = Vis(vis_file)
plt=nothing
#=
vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end
=#

include("../examples/remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

#manually set the initial yaw angles for the optimisation 
#con.yaw_data = [sim.start_time:sim.end_time        180 .*  ones(sim.end_time-sim.start_time+1, 9)]


# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();


wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);

#run initial conditions
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 




#------------------------------------------------------------------------------------------------------------
#define required cost and memory helper functions


function construct_yaw_matrix(x, sim, wf)   #wf to be used later


    yaws =  ones(sim.end_time - sim.start_time + 1)   * x' * 360.0  #expands into a matrix,
                                                                    # x is in [0,1] range
        
    #ADD CHECK: yaws must not excees X degrees misalignment from wind to prevent crash
    return   [sim.start_time:sim.end_time    yaws]
end




function apply_yaw_rate_limit!(yaws, max_yaw_rate)
    ax1, ax2 = axes(yaws)
    
    @inbounds for i in Iterators.drop(ax1, 1)
        for j in ax2
            yaw_change = yaws[i, j] - yaws[i-1, j]
            
            if abs(yaw_change) > max_yaw_rate
                num_steps = ceil(Int, abs(yaw_change) / max_yaw_rate)
                yaw_step = yaw_change / num_steps
                
                for k in 1:num_steps
                    yaws[i-k+1, j] = yaws[i-k, j] + yaw_step
                end
            end
        end
    end
end



#dynamic yaw matrix construction, where the optimiser can choose at which time steps the yaw angles change, and what the yaw angles are for each time period between changes
function construct_yaw_matrix_dynamic(x, sim, wf, num_yaw_changes, max_yaw_rate)   #wf to be used later

    #structure of x: [time_yaw_change | yaw_change_vector1 | yaw_change_vector2 | ...]
    # where time_yaw_change is a vector of length num_yaw_changes-1, with values
    #I want to to let the optimiser choose at which time steps the yaw angles change, so these are the time_yaw_change variables.
    # the yaw_change_vectors are the actual yaw angles for each change, which are held constant between the time steps defined in time_yaw_change
    # so, for example, if num_yaw_changes = 3, then x would look like this: [t1, t2, yaw1_vector, yaw2_vector, yaw3_vector]
    # where t1 and t2 are the time steps at which the yaw angles change, and yaw1_vector, yaw2_vector, yaw3_vector are the yaw angles for each change, which are held constant between the time steps defined in time_yaw_change
    

    if num_yaw_changes == 1
        yaws = ones(sim.end_time - sim.start_time + 1)   * x' * 360.0  #expands into a matrix,
                                                                    # x is in [0,1] range
    
        return  [sim.start_time:sim.end_time    yaws]
        #return statement ends function
    end


    #do this preallocation earlier and then fill it, probably also inplace writing
    yaws = Matrix{Float64}(undef, length(sim.start_time:sim.end_time), wf.nT)  #initialize empty matrix to store yaw angles for each time step



    yaw_change_timestamps = [0   x[1:num_yaw_changes-1]'] .* (sim.end_time - sim.start_time) .+ sim.start_time   #scale back to time range
    yaw_change_timestamps = round.(Int, yaw_change_timestamps)  #round to integers for indexing
    yaw_changes = reshape(x[num_yaw_changes:end], (wf.nT, num_yaw_changes))' 
    #reshape yaw changes into matrix of size [num_yaw_changes x nT]
    yaws = ones(yaw_change_timestamps[2] - yaw_change_timestamps[1] + 1)  * yaw_changes[1, :]' .* 360.0  #first time period, from start to first change
    #                                                        +1 for the 0th data entry
    
    
    for i = 2:num_yaw_changes-1
        #set yaws of current step
        yaw =  ones(yaw_change_timestamps[i+1] - yaw_change_timestamps[i])  *  yaw_changes[i, :]' .* 360.0  #expands into a matrix, where each row is the yaw angles for that time period
        yaws = [yaws; yaw]
    end

    #last time period, from last change to end
    yaw =  ones(sim.end_time - yaw_change_timestamps[end])  *  yaw_changes[end, :]' .* 360.0 

    yaws = [yaws; yaw]
    #note: fix allocation for memory bottleneck
    

    # limit transition to max_yaw_rate, by linearly interpolating between the yaw angles before and after the change over a time period defined by max_yaw_rate

    apply_yaw_rate_limit!(yaws, max_yaw_rate)  #in-place modification of yaw matrix to limit yaw rate,
                                                #by linearly interpolating between the yaw angles
                                                #before and after the change over a time 
                                                #period defined by max_yaw_rate

    #plot yaw over time ...
    return   [sim.start_time:sim.end_time    yaws ]
end



function create_fitness(plt, set, wf::WindFarm, wind::Wind, sim, con, vis, floridyn, floris)
    # Create task-local copies of mutable state objects
    # These get deep-copied once per task, ensuring no shared state
    tlv = TaskLocalValue{NamedTuple}() do
        (
            
            plt = deepcopy(plt),
            set = deepcopy(set),
            wf = deepcopy(wf),
            wind = deepcopy(wind),
            sim = deepcopy(sim),
            floris = deepcopy(floris),
            floridyn = deepcopy(floridyn),
            vis = deepcopy(vis),
            con = deepcopy(con)
            #need to figure out which input arguments get modified
            # in-place to determine what needs to be deep-copied
        )
    end
    
    return function cost(x)
 
          
        
        
        state = tlv[]
        # Call runFLORIDyn with task-local copies
        # So each thread modifies its own copies, never shared
        plt_local  = state.plt
        set_local  = state.set
        wf_local = state.wf
        wind_local = state.wind
        sim_local = state.sim
        floris_local = state.floris
        floridyn_local = state.floridyn
        vis_local = state.vis
        con_local = state.con

        #mapping from x to yaw matrix
    
        con_local.yaw_data = construct_yaw_matrix_dynamic(x, sim_local, wf_local, set_num_yaw_changes, set_max_yaw_rate)  #construct yaw matrix for current candidate solution x, with max yaw rate of 10 deg/s



        # Now safe: each thread has isolated state
       
        wf, md, mi = run_floridyn(plt_local, set_local, wf_local, wind_local, sim_local, con_local, vis_local, floridyn_local, floris_local)
        return -sum(md.PowerGen)
        

        #could try something like a try/catch statement to expand search space, 
        #and gradient towards wind direction to return simulation back to feasible direction
        #instead of providing a flat indicator function

        #...
        

        
    end
end




"""
 function construct_yaw_matrix()
        #construct yaw matrix for optimisation, reducing the number of free variables for the optimiser


"""



#------------------------------------------------------------------------------------------------------------
#OPTIMISATION PART


#number of yaw changes allowed
set_num_yaw_changes = 3 #N
set_max_yaw_rate = 1.0 #deg/s

# set cost function
fit_func = create_fitness(plt, set, wf, wind, sim, con, vis, floridyn, floris)


#init: yaws aligned with wind, equal time spacing
equal_time_spacing = 0:1/set_num_yaw_changes:1
#include wind direction reading here
yaw_guesses = repeat(0/360 * ones(wf.nT), set_num_yaw_changes)  #do some interpolation on 
                                # wind direction here






#set initial guess

#x0 = vcat(equal_time_spacing[2:end-1], yaw_guesses)  #182 deg
x0 = result.minimizer  #start from previous result








#AAAAAAA note for next time: make these limits dependent on the dynamic wind field
lower_bounds = vcat(zeros(set_num_yaw_changes-1),    -0.4 * ones(wf.nT))
upper_bounds = vcat(ones(set_num_yaw_changes-1),      0.4 * ones(wf.nT))






opts = Evolutionary.Options(
    iterations = 30,
    abstol = 1e-8,
    reltol = 1e-8,
    show_trace = true,
    show_every = 1,
    store_trace = true,
    parallelization = :thread  #serial   #thread   
)


#max yawing rate



#hyperparams
set_lambda_multiplier = 3
set_lambda0 = 2 * round(   (4 + 3 * log(wf.nT)) * set_num_yaw_changes     / 2.0   )   #half the offspring 
set_lambda = Int(set_lambda_multiplier * set_lambda0)
set_mu = Int(round(set_lambda / 2))
set_sigma0 =  0.03  # set to 30% of the search range


#think about rescaling for time steps, this sigma is too small that part of the optim space



@time result =  Evolutionary.optimize(fit_func,
                                BoxConstraints(lower_bounds, upper_bounds),
                                x0, 
                                CMAES(  lambda = set_lambda,
                                        mu = set_mu,
                                        sigma0 = set_sigma0),
                                opts)






#----------------------------------------------------------------------------------

#PLOTTING







# the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)


#set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set = Settings(wind, sim, con, false, false)
set.enable_veer = true
set.control_mode = Yaw_Optimisation();

#run initial conditions
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim);
wf = initSimulation(wf, sim);

#disable online visualisation
vis.online = false 


#plot flowfield
con.yaw_data = construct_yaw_matrix_dynamic(result.minimizer, sim, wf, set_num_yaw_changes, set_max_yaw_rate)
wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)


#include("get_energy.jl")  #for energy tests
#include("controller_output.jl")   #to see and plot controller angles






#=
AAAAAAA        
#the idea is now to limit the yaw angles explored by the optimiser: limit them from max.
#+/- 90 deg from current wind direction








I think something is wrong with the multithreading is a bit off,
as the globally best result is not returned while calling get_energy.jl

#include automatic IPOP restart strategies for CMA-ES to get out of local minima?






=#