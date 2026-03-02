#the include file for all functions


#------------------------------------------------------------------------------------------------------------
#define required cost and memory helper functions

#old implementation for static yaw control
function construct_yaw_matrix(x, sim, wf)   #wf to be used later


    yaws =  ones(sim.end_time - sim.start_time + 1)   * x' * 360.0  #expands into a matrix,
                                                                    # x is in [0,1] range
        
    #ADD CHECK: yaws must not excees X degrees misalignment from wind to prevent crash
    return   [sim.start_time:sim.end_time    yaws]
end


#time-dependent yaw matrix construction
function construct_yaw_matrix_dynamic(x, sim, wf, num_yaw_changes, max_yaw_rate)
    
    num_steps = sim.end_time - sim.start_time + 1
    
    # 1. Preallocate the final yaw matrix exactly once to eliminate memory bottlenecks
    yaws = Matrix{Float64}(undef, num_steps, wf.nT)
    
    if num_yaw_changes == 1
        # Fill directly using broadcasting instead of allocating 'ones' matrices
        for j in 1:wf.nT
            yaws[:, j] .= x[j] * 360.0
        end
        return hcat(sim.start_time:sim.end_time, yaws)
    end

    duration = sim.end_time - sim.start_time
    
    # 2. Map normalized timestamps to specific row boundaries in our preallocated matrix
    transitions = zeros(Int, num_yaw_changes)
    @inbounds for i in 1:(num_yaw_changes - 1)
        # Calculate row index corresponding to the time change
        transitions[i] = round(Int, x[i] * duration) + 1
    end
    transitions[end] = num_steps # Ensure the last segment reaches the very end
    
    # 3. Fill the matrix in-place
    current_row = 1
    for i in 1:num_yaw_changes
        end_row = transitions[i]
        
        # Safeguard: optimizers sometimes guess non-sequential times (e.g., t2 < t1)
        end_row = max(current_row - 1, end_row) 
        if i == num_yaw_changes
            end_row = num_steps # Force last segment to cover the remaining steps
        end
        
        if current_row <= end_row
            @inbounds for j in 1:wf.nT
                # Read directly from flat array 'x' using an index formula.
                # This entirely avoids expensive `reshape` and array slicing.
                flat_idx = num_yaw_changes - 1 + (i - 1) * wf.nT + j
                val = x[flat_idx] * 360.0
                
                # Assign the scalar value directly to the specific block of memory
                yaws[current_row:end_row, j] .= val 
            end
        end
        current_row = end_row + 1
    end

    # Enforce yaw rate limits (in-place modification)
    apply_yaw_rate_limit!(yaws, max_yaw_rate)

    # hcat binds the time vector and the yaw matrix without generating extra rows
    return hcat(sim.start_time:sim.end_time, yaws)
end


#parallel test function without try/catch to see errors
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

#yaw rate limiting
function apply_yaw_rate_limit!(yaws, max_yaw_rate)
    ax1, ax2 = axes(yaws)
    
    @inbounds for i in Iterators.drop(ax1, 1)
        @inbounds for j in ax2
            yaw_change = yaws[i, j] - yaws[i-1, j]
            
            if abs(yaw_change) > max_yaw_rate
         
                num_steps = ceil(Int, abs(yaw_change) / max_yaw_rate)
                yaw_step = yaw_change / num_steps
                
                @inbounds for k in 1:num_steps
                    yaws[i-k+1, j] = yaws[i-k, j] + yaw_step
                end
            end
        end
    end
end

#wind direction interpolation helper function
function get_wind_at_t(t, wind_matrix)
    # wind_matrix is [time direction]
    times = wind_matrix[:, 1]
    dirs  = wind_matrix[:, 2]
    
    # Handle bounds
    if t <= times[1];   return dirs[1];   end
    if t >= times[end]; return dirs[end]; end
    
    # Find the interval for linear interpolation
    idx = findfirst(x -> x >= t, times)
    t_low, t_high = times[idx-1], times[idx]
    d_low, d_high = dirs[idx-1], dirs[idx]
    
    # Linear interpolation formula
    return d_low + (d_high - d_low) * (t - t_low) / (t_high - t_low)
end

#generates an initial guess aligned with the middle point of the uniform time segments
function generate_initial_guess(sim, wind, wf, n_segments)
    duration = sim.end_time - sim.start_time
    # equal_time_spacing gives the boundaries: [t0, t1, t2... tn]
    equal_time_spacing = collect(0:1/n_segments:1)
    
    # 1. Internal transition times (normalized 0 to 1)
    # These are the 'knots' the optimizer moves to change segment lengths
    time_guesses = equal_time_spacing[2:end-1]
    
    # 2. Yaw guesses (normalized degrees / 360)
    yaw_guesses = Float64[]
    
    for i in 1:n_segments
        # Find the window boundaries in normalized time
        t_start_norm = equal_time_spacing[i]
        t_end_norm   = equal_time_spacing[i+1]
        
        # Calculate the midpoint of the window to get the "average" interpolated wind
        t_mid_norm = (t_start_norm + t_end_norm) / 2.0
        
        # Scale to actual simulation time
        t_actual = duration * t_mid_norm + sim.start_time
        
        # Get interpolated wind at the midpoint
        wind_dir = get_wind_at_t(t_actual, wind.dir)
        
        # Normalize for the optimizer (x * 360 = degrees)
        val_norm = wind_dir / 360.0
        
        # Fill the guess for all turbines in this window
        append!(yaw_guesses, fill(val_norm, wf.nT))
    end
    
    return vcat(time_guesses, yaw_guesses)
end


@inline function l1_norm_penalty(yaws)
    rows, cols = size(yaws)
    total_sum = 0.0
    
    # Iterate column-first (Julia is column-major, this is faster)
    @inbounds for j in 1:cols
        for i in 1:(rows - 1)
            # Accessing indices directly avoids creating slices or intermediate arrays
            total_sum += abs(yaws[i, j] - yaws[i+1, j])
        end
    end
    
    return total_sum / wf.nT / (sim.end_time - sim.start_time + 1)
end



#cost function structure with parallel, and constraint handling
function parallel_costfunction(plt, set, wf::WindFarm, wind::Wind, sim, con, vis, floridyn, floris, objective)
    
    tlv = TaskLocalValue{NamedTuple}() do
        (
            plt=deepcopy(plt), set=deepcopy(set), wf=deepcopy(wf),
            wind=deepcopy(wind), sim=deepcopy(sim), floris=deepcopy(floris),
            floridyn=deepcopy(floridyn), vis=deepcopy(vis), con=deepcopy(con),
            general_purpose_buffer = zeros(wf.nT)
        )
    end
    
    return function cost(x)
        state = tlv[]
        penalty_term = 0.0
    

        #yaw matrix construction
        
        yaws_with_time = construct_yaw_matrix_dynamic(x, state.sim, state.wf, set_num_yaw_changes, set_max_yaw_rate)
        
        #get data for penalty function
        sim_times = @view yaws_with_time[:, 1]
        yaws_only = @view yaws_with_time[:, 2:end] # All turbines
        wind_table = state.wind.dir          
        
        #check yaw alignment
        max_violation = 0.0
        @inbounds for (i, t) in enumerate(sim_times)
            current_wind_dir = get_wind_at_t(t, wind_table)
          
            # Check all turbines at this specific second
            @inbounds for turbine_idx in 1:wf.nT
                rel_yaw = abs(yaws_only[i, turbine_idx] - current_wind_dir)
                # Keep track of the worst offender across all turbines and time steps
                max_violation = max(max_violation, rel_yaw - set_max_yaw_misalignment)
            end
        end

        #check maximum yaw l1 of change, of any turbine
        # Check maximum yaw L1 change (Zero Allocations)
        num_rows = size(yaws_only, 1)
        
        @inbounds for j in 1:wf.nT # For each turbine (column)
            acc = 0.0
            @inbounds for i in 1:(num_rows - 1)
                # Access indices directly to avoid slicing and temporary arrays
                acc += abs(yaws_only[i, j] - yaws_only[i+1, j])
            end
            
            # Store the result in our preallocated buffer from TaskLocalValue
            state.general_purpose_buffer[j] = acc
            
            if acc > set_lambda_l1_hard_limit 
                return 1e6 + (acc - set_lambda_l1_hard_limit)^2 * 1000.0
            end
        end



        #in case turbines are misaligned too strongly  
        if max_violation > 0
            # Quadratic penalty creates a smooth "slope" leading back to 0 violation
            return 1e6 + (max_violation^2 * 1000.0)
        end

        #any other constraint code:
        penalty_term = penalty_term + set_lambda_l1 * l1_norm_penalty(yaws_only)



        #otherwise
        try
            state.con.yaw_data = yaws_with_time
            wf, md, mi = run_floridyn(
                state.plt, state.set, state.wf, state.wind, 
                state.sim, state.con, state.vis, state.floridyn, state.floris
            )
            
            
            #objective

            return objective(md) + penalty_term  #in kW per turbine
            #return powerTrackingObjective(md) + penalty_term



        catch e
            # 1. If the user hits Ctrl+C, let it happen!
            if e isa InterruptException 
                rethrow(e)
            end
            if e isa ArgumentError #memory errors for Evolutionary.jl, fall under this
                rethrow(e)
            end
            return 2e6 # Fallback for unexpected failures
        end
    end
end






#objectivefunctions
@inline function totalEnergyObjective(md)

    return -sum(md.PowerGen) / (wf.nT * sim.n_sim_steps) * 1000.0   #in kW per turbine
    
end


function powerTrackingObjective(md)
    #some kind of way to track a power curve by returning the l2 norm of the difference between the power curve and the actual power generated, normalized by the number of turbines and time steps
    #(and plot... add plotting scripts) 

    power_per_step = zeros(sim.n_sim_steps)
    for k in 1:sim.n_sim_steps
    acc = 0.0
    offset = (k - 1) * wf.nT
    @inbounds for j in 1:wf.nT
        acc += md.PowerGen[offset + j]
            end
    power_per_step[k] = acc
    end
    
    return norm(power_per_step - set_desired_power_curve) / (wf.nT * sim.n_sim_steps) * 1000 #rms error in kW per turbine
end








"""
this set of functions allows for individual turbine switching times,
but this makes the optim. problem too hard, and having a common switching time
leads to better results.

"""



# cost function with for individual turbine switching plus trycatch for misalignment
function parallel_costfunction_individual_switches(plt, set, wf::WindFarm, wind::Wind, sim, con, vis, floridyn, floris, objective)
    tlv = TaskLocalValue{NamedTuple}() do
        (
            plt=deepcopy(plt), set=deepcopy(set), wf=deepcopy(wf),
            wind=deepcopy(wind), sim=deepcopy(sim), floris=deepcopy(floris),
            floridyn=deepcopy(floridyn), vis=deepcopy(vis), con=deepcopy(con),
            general_purpose_buffer = zeros(wf.nT)
        )
    end
    
    return function cost(x)
        state = tlv[]
        penalty_term = 0.0
    

        #yaw matrix construction
        
        yaws_with_time = construct_yaw_matrix_dynamic_individual_switches(x, state.sim, state.wf, set_num_yaw_changes, set_max_yaw_rate)
        
        #get data for penalty function
        sim_times = @view yaws_with_time[:, 1]
        yaws_only = @view yaws_with_time[:, 2:end] # All turbines
        wind_table = state.wind.dir          
    
        #check yaw alignment
        max_violation = 0.0
        @inbounds for (i, t) in enumerate(sim_timeAAAAAAAA)  
            #AAAAAAAAAAAAAAAA    note: crashing caused by wrong time retrieval
            #detect changes from delta yaw matrix
            #retrieve timestamps from there
            #then insert those timestamps into the code below to detects yaw misalignment more accurately
            current_wind_dir = get_wind_at_t(t, wind_table)
          
            # Check all turbines at this specific second
            @inbounds for turbine_idx in 1:wf.nT
                rel_yaw = abs(yaws_only[i, turbine_idx] - current_wind_dir)
                # Keep track of the worst offender across all turbines and time steps
                max_violation = max(max_violation, rel_yaw - set_max_yaw_misalignment)
            end
        end

        #check maximum yaw l1 of change, of any turbine
        # Check maximum yaw L1 change (Zero Allocations)
        num_rows = size(yaws_only, 1)
        
        @inbounds for j in 1:wf.nT # For each turbine (column)
            acc = 0.0
            @inbounds for i in 1:(num_rows - 1)
                # Access indices directly to avoid slicing and temporary arrays
                acc += abs(yaws_only[i, j] - yaws_only[i+1, j])
            end
            
            # Store the result in our preallocated buffer from TaskLocalValue
            state.general_purpose_buffer[j] = acc
            
            if acc > set_lambda_l1_hard_limit 
                return 1e6 + (acc - set_lambda_l1_hard_limit)^2 * 1000.0 #make sure this is never lower than actual objective
            end
        end



        #in case turbines are misaligned too strongly  
        if max_violation > 0
            # Quadratic penalty creates a smooth "slope" leading back to 0 violation
            return 1e6 + (max_violation^2 * 1000.0)
        end

        #any other constraint code:
        penalty_term = penalty_term + set_lambda_l1 * l1_norm_penalty(yaws_only)



        #otherwise
        try
            state.con.yaw_data = yaws_with_time
            wf, md, mi = run_floridyn(
                state.plt, state.set, state.wf, state.wind, 
                state.sim, state.con, state.vis, state.floridyn, state.floris
            )
            
            
            #objective

            return objective(md) + penalty_term  #in kW per turbine
            #return powerTrackingObjective(md) + penalty_term



        catch e
            # 1. If the user hits Ctrl+C, let it happen!
            if e isa InterruptException 
                rethrow(e)
            end
            if e isa ArgumentError #memory errors for Evolutionary.jl, fall under this
                rethrow(e)
            end
            return 2e6 # Fallback for unexpected failures
        end
    end
end



# time-dependent yaw matrix construction (Individual Turbine Switching)
function construct_yaw_matrix_dynamic_individual_switches(x, sim, wf, num_yaw_changes, max_yaw_rate)
    
    num_steps = Int(sim.end_time - sim.start_time + 1)
    duration = sim.end_time - sim.start_time
    num_ts_params = num_yaw_changes - 1
    
    # Preallocate the final yaw matrix
    yaws = Matrix{Float64}(undef, num_steps, wf.nT)
    
    # Calculate the global offset where the yaw values begin
    # Total time entries = (number of changes - 1) * number of turbines
    yaw_start_offset = num_ts_params * wf.nT

    for j in 1:wf.nT
        # 1. Extract transition times for turbine j
        # These are located at [(j-1)*num_ts_params + 1] to [j*num_ts_params]
        t_indices = zeros(Int, num_yaw_changes)
        time_block_start = (j - 1) * num_ts_params
        
        for i in 1:num_ts_params
            t_val = x[time_block_start + i]
            t_indices[i] = round(Int, t_val * duration) + 1
        end
        t_indices[end] = num_steps # Final segment boundary

        # 2. Fill the column for turbine j
        current_row = 1
        yaw_block_start = yaw_start_offset + (j - 1) * num_yaw_changes
        
        for i in 1:num_yaw_changes
            # Extract yaw value from the second half of vector x
            val = x[yaw_block_start + i] * 360.0
            
            end_row = t_indices[i]
            
            # Robustness: ensure segments don't overlap backward
            end_row = clamp(end_row, current_row - 1, num_steps)
            
            if current_row <= end_row
                @inbounds yaws[current_row:end_row, j] .= val
            end
            current_row = end_row + 1
        end
    end

    # Enforce yaw rate limits (in-place modification)
    apply_yaw_rate_limit!(yaws, max_yaw_rate)
    return hcat(sim.start_time:sim.end_time, yaws)

end


# Generates an initial guess with individual turbine switching times
function generate_initial_guess_individual_switches(sim, wind, wf, n_segments)
    duration = sim.end_time - sim.start_time
    nT = wf.nT
    
    # equal_time_spacing gives the boundaries: [t0, t1, t2... tn]
    equal_time_spacing = collect(0:1/n_segments:1)
    
    # 1. Individual transition times (normalized 0 to 1)
    # We now need (n_segments - 1) * nT time entries
    time_guesses = Float64[]
    
    # We start from the second element of equal_time_spacing (the first transition)
    # and go up to the second to last (the last internal transition)
    for i in 2:n_segments
        t_norm = equal_time_spacing[i]
        # Repeat the same transition time for every turbine as the starting point
        append!(time_guesses, fill(t_norm, nT))
    end
    
    # 2. Yaw guesses (normalized degrees / 360)
    # We need n_segments * nT yaw entries
    yaw_guesses = Float64[]
    
    for i in 1:n_segments
        # Calculate the midpoint of the window
        t_mid_norm = (equal_time_spacing[i] + equal_time_spacing[i+1]) / 2.0
        t_actual = duration * t_mid_norm + sim.start_time
        
        # Get interpolated wind at the midpoint
        wind_dir = get_wind_at_t(t_actual, wind.dir)
        val_norm = wind_dir / 360.0
        
        # Fill the guess for all turbines in this segment
        append!(yaw_guesses, fill(val_norm, nT))
    end
    
    # Combine into the final flat vector x
    return vcat(time_guesses, yaw_guesses)
end



#old: --------------------------------------------------------------------------------------



#"testing area inside testing area"----------------------------------------------------------------------------


@inline function l1_norm_penalty_slow(yaws)
    # Calculate the L1 norm of the yaw angles, currently farm total
    return sum(abs.(yaws[1:end-1, :] - yaws[2:end, :] )) / wf.nT / (sim.end_time - sim.start_time + 1)  # Normalize by number of turbines and time
end










