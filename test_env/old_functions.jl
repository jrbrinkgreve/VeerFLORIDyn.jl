
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




#------------------------------------------------------------------------------------------

#19 feb 16:00 before adding trycatch statement backup




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




#3 march : replaced function by allocationless version
#time-dependent yaw matrix construction



#------------------------------------------------------------------------------------------

#19 feb 16:40 replacing with the automatix wind direction aligned initialisation function





#=
#init: yaws aligned with wind, equal time spacing
equal_time_spacing = 0:1/set_num_yaw_changes:1
#include wind direction reading here

yaw_guesses = repeat(0/360 * ones(wf.nT), set_num_yaw_changes)  #do some interpolation on 
                                # wind direction here
#set initial guess
x0 = vcat(equal_time_spacing[2:end-1], yaw_guesses)  #time; deg
#x0 = result.minimizer  #start from previous result
=#




