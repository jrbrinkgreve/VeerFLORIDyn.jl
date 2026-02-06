using Infiltrator

function construct_yaw_matrix(x, sim, wf, num_yaw_changes)   #wf to be used later

    #structure of x: [time_yaw_change | yaw_change_vector1 | yaw_change_vector2 | ...]
    # where time_yaw_change is a vector of length num_yaw_changes-1, with values
    #I want to to let the optimiser choose at which time steps the yaw angles change, so these are the time_yaw_change variables.
    # the yaw_change_vectors are the actual yaw angles for each change, which are held constant between the time steps defined in time_yaw_change
    # so, for example, if num_yaw_changes = 3, then x would look like this: [t1, t2, yaw1_vector, yaw2_vector, yaw3_vector]
    # where t1 and t2 are the time steps at which the yaw angles change, and yaw1_vector, yaw2_vector, yaw3_vector are the yaw angles for each change, which are held constant between the time steps defined in time_yaw_change
    
    
    yaws = Matrix{Float64}(undef, length(sim.start_time:sim.end_time), wf.nT)  #initialize empty matrix to store yaw angles for each time step

    yaw_change_timestamps = [0   x[1:num_yaw_changes-1]'] .* (sim.end_time - sim.start_time) .+ sim.start_time   #scale back to time range
    yaw_change_timestamps = round.(Int, yaw_change_timestamps)  #round to integers for indexing
    yaw_changes = reshape(x[num_yaw_changes:end], (wf.nT, num_yaw_changes))' 
    #reshape yaw changes into matrix of size [num_yaw_changes x nT]
    yaws = ones(yaw_change_timestamps[2] - yaw_change_timestamps[1] + 1)  * yaw_changes[1, :]' .* 360.0  #first time period, from start to first change
    #                                                        +1 for the 0th data entry
    
    
    for i = 2:num_yaw_changes-1


        #note: possible error with timestamps variable creating arrays of different sizes
        yaw =  ones(yaw_change_timestamps[i+1] - yaw_change_timestamps[i])  *  yaw_changes[i, :]' .* 360.0  #expands into a matrix, where each row is the yaw angles for that time period
        yaws = [yaws; yaw]
    end
    #last time period, from last change to end
    yaw =  ones(sim.end_time - yaw_change_timestamps[end])  *  yaw_changes[end, :]' .* 360.0

    yaws = [yaws; yaw]

    #ADD CHECK: yaws must not excees X degrees misalignment from wind to prevent crash
    return   [sim.start_time:sim.end_time    yaws]
end




#=
example:

x = [0.3 0.6 | 0.5 0.4 0.5 0.4 | 0.2 0.3 0.2 0.3 | 0.1 0.1 0.1 0.1]

the first 2 values are the time_yaw_change variables, which define the time steps at which the yaw angles change.
So in this example, the yaw angles would change at 30% and 60% of the simulation time.

the next 3 blocks of 4 values are the yaw_change_vectors, which define the yaw angles for each change, which are held constant between the time steps defined in time_yaw_change.
So in this example, the yaw angles would be 0.5, 0.4, 0.5, 0.4 for the first time period
(from start to 30% of the simulation time), then 0.2, 0.3, 0.2, 0.3 for the second time period
(from 30% to 60% of the simulation time), and then 0.1, 0.1, 0.1, 0.1 for the third time period
(from 60% to the end of the simulation time).


=#


#construct_yaw_matrix([0.25 0.75 0.5 0.4 0.5 0.4  0.2 0.3 0.2 0.3 0.1 0.1 0.1 0.1], sim, wf, 3)




for i = 0:1/10:1
    println(i)
end
