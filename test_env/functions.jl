#the include file for all functions


#------------------------------------------------------------------------------------------------------------
#define required cost and memory helper functions



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


function construct_yaw_matrix_dynamic!(buffer, x, sim, wf, opt_set)
    num_yaw_changes = opt_set.set_num_yaw_changes
    max_yaw_rate = opt_set.set_max_yaw_rate
    num_steps = size(buffer, 1)
    duration = sim.end_time - sim.start_time
    
    # 1. Fill the first column with timestamps (in-place)
    @inbounds for j in 1:num_steps
        buffer[j, 1] = sim.start_time + (j - 1)
    end

    # Handle the simple case (1 yaw change = constant yaw)
    if num_yaw_changes == 1
        for j in 1:wf.nT
            val = x[j] * 360.0
            col_idx = j + 1
            @inbounds for i in 1:num_steps
                buffer[i, col_idx] = val
            end
        end
        # Enforce rate limits on the yaw portion (columns 2 to end)
        apply_yaw_rate_limit!(@view(buffer[:, 2:end]), max_yaw_rate)
        return nothing
    end

    # 2. Process segments
    current_row = 1
    for i in 1:num_yaw_changes
        # Determine end_row for this time segment without allocating a 'transitions' vector
        if i < num_yaw_changes
            # x[i] is the normalized timestamp
            end_row = round(Int, x[i] * duration) + 1
        else
            end_row = num_steps
        end
        
        # Safeguard: Ensure end_row doesn't go backwards
        end_row = clamp(end_row, current_row, num_steps)
        
        if current_row <= end_row
            for j in 1:wf.nT
                # Calculate the flat index for the yaw value in x
                # (Num_time_params) + (segment_idx * nT) + turbine_idx
                flat_idx = (num_yaw_changes - 1) + (i - 1) * wf.nT + j
                val = x[flat_idx] * 360.0
                
                col_idx = j + 1
                @inbounds for row_idx in current_row:end_row
                    buffer[row_idx, col_idx] = val
                end
            end
        end
        current_row = end_row + 1
        if current_row > num_steps break end
    end

    # 3. Enforce yaw rate limits in-place on the yaw columns
    yaws_view = @view buffer[:, 2:end]
    apply_yaw_rate_limit!(yaws_view, max_yaw_rate)
    
    return nothing
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



function fill_wind_dir_buffer!(buffer, sim_times, wind_matrix)
    # dest_buffer: 1201 elements
    # sim_times: 1201 elements
    # wind_matrix: [Time Direction]
    
    w_times = @view wind_matrix[:, 1]
    w_dirs  = @view wind_matrix[:, 2]
    n_wind  = length(w_times)
    curr_idx = 2 # Start at the first possible high interval
    
    @inbounds for i in eachindex(sim_times)
        t = sim_times[i]
        
        # Clamp bounds
        if t <= w_times[1]
            buffer[i] = w_dirs[1]
            continue
        elseif t >= w_times[end]
            buffer[i] = w_dirs[end]
            continue
        end
        
        # Advance the wind_matrix index only as far as needed (Sweeping algorithm)
        while curr_idx < n_wind && w_times[curr_idx] < t
            curr_idx += 1
        end
        
        # Interpolate
        t_low,  t_high = w_times[curr_idx-1], w_times[curr_idx]
        d_low,  d_high = w_dirs[curr_idx-1],  w_dirs[curr_idx]
        
        buffer[i] = d_low + (d_high - d_low) * (t - t_low) / (t_high - t_low)
    end
    return nothing
end


#main cost function
function parallel_costfunction(plt, set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, vis::Vis, floridyn::FloriDyn, floris::Floris, opt_set::OptimisationSettings)
    
    return function cost(x)
        state = tlv[]
        penalty_term = 0.0
    

        #yaw matrix construction
        construct_yaw_matrix_dynamic!(state.con.yaw_data, x, state.sim, state.wf, opt_set)
        
        #get data for penalty function
        sim_times = @view state.con.yaw_data[:, 1]
        yaws_only = @view state.con.yaw_data[:, 2:end] # All turbines   
        
        #get full wind direction matrix
        fill_wind_dir_buffer!(state.wind_dirs_buffer, sim_times, state.wind.dir)


        #handle soft and hard constraints
        penalty_term = calculate_penalties(yaws_only, state.wind_dirs_buffer, state.sim, state.wf.nT, opt_set)
        
        #detect if penalty term exceeds hard limit, and if so return it directly to avoid unnecessary simulations
        if penalty_term >= 1e6
            return penalty_term
        end
        
        #otherwise
        try
            runFLORIDyn_optimisation!(state.power_vector, state.plt, state.set, state.wf, state.wind, 
                state.sim, state.con, state.vis, state.floridyn, state.floris)
            
            #objective
            return totalEnergyObjective(state.power_vector, state.wf.nT, state.sim.n_sim_steps) + penalty_term  #in kW per turbine

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


#constrains handling
function calculate_penalties(yaws_only, wind_dirs, sim::Sim, nT::Int, opt_set::OptimisationSettings)
    max_violation = 0.0
    total_l1_acc  = 0.0 

    # Column-Major traversal (Turbines then Time) for maximum cache efficiency
    @inbounds for i in 1:nT             
        turbine_l1 = 0.0
        
        @inbounds for j in 1:(sim.end_time - sim.start_time + 1)       
            curr_yaw = yaws_only[j, i]
            
            # 1. Check Alignment Violation
            v = abs(curr_yaw - wind_dirs[j]) - opt_set.set_max_yaw_misalignment
            max_violation = max(max_violation, v)

            # 2. Accumulate L1 change
            if j > 1
                prev_yaw = yaws_only[j-1, i]
                turbine_l1 += abs(curr_yaw - prev_yaw)
            end
        end
        
        # 3. Early Exit for L1 Hard Limit
        if turbine_l1 > opt_set.set_lambda_l1_hard_limit
            return 1e6 + (turbine_l1 - opt_set.set_lambda_l1_hard_limit)^2 * 1000.0
        end
        
        total_l1_acc += turbine_l1
    end

    # 4. Early Exit for Misalignment Hard Limit
    if max_violation > 0
        return 1e6 + (max_violation^2 * 1000.0)
    end

    # 5. Return normalized soft penalty
    return (opt_set.set_lambda_l1 * total_l1_acc) / (nT * (sim.end_time - sim.start_time + 1))
end



#objective functions
@inline function totalEnergyObjective(powervector, nT, nsteps)

    return -sum(powervector) / (nT * nsteps) * 1000.0   #in kW per turbine
    
end

#wind direction interpolation helper function
#replaced by fill_wind_dir_buffer! for performance
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








#unused below:



#power tracking ---------------------------------------------------------------------------------
function powerTrackingObjective(md, sim, wf, set_desired_power_curve)
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


# cost function with for individual turbine switching plus trycatch for misalignment
function parallel_costfunction_individual_switches(plt, set, wf::WindFarm, wind::Wind, sim, con, vis, floridyn, floris, objective)
    tlv = TaskLocalValue{NamedTuple}() do
        (
            plt=deepcopy(plt), set=deepcopy(set), wf=deepcopy(wf),
            wind=deepcopy(wind), sim=deepcopy(sim), floris=deepcopy(floris),
            floridyn=deepcopy(floridyn), vis=deepcopy(vis), con=deepcopy(con),

            wind_dirs_buffer = zeros(sim.end_time - sim.start_time + 1), #preallocate buffer for wind directions at each time step to avoid allocations in get_wind_at_t
            yaws_with_time_buffer = zeros(sim.end_time - sim.start_time + 1, wf.nT + 1) #preallocate buffer for yaw matrix to avoid allocations in construct_yaw_matrix_dynamic_individual_switches
        )
    end
    
    return function cost(x)
        state = tlv[]
        penalty_term = 0.0
    

        #yaw matrix construction
        construct_yaw_matrix_dynamic_individual_switches!(state.yaws_with_time_buffer, x, state.sim, state.wf, set_num_yaw_changes, set_max_yaw_rate)
        
        #get data for penalty function
        sim_times = @view state.yaws_with_time_buffer[:, 1]
        yaws_only = @view state.yaws_with_time_buffer[:, 2:end] # All turbines
                
    
        state.wind_dirs_buffer .= get_wind_at_t.(sim_times, Ref(state.wind.dir))        

        max_violation = 0.0
        total_l1_acc  = 0.0

        # Performance: Loop over Columns (Turbines) then Rows (Time)
        # This accesses memory sequentially (Column-Major)
        @inbounds for i in 1:wf.nT             
            turbine_l1 = 0.0
            
            for j in 1:(sim.end_time - sim.start_time + 1)       
                curr_yaw = yaws_only[j, i]
                
                # Check Alignment Violation relative to wind at time j
                v = abs(curr_yaw - state.wind_dirs_buffer[j]) - set_max_yaw_misalignment
                max_violation = max(max_violation, v)

                # Accumulate L1 change for this turbine
                if j > 1
                    prev_yaw = yaws_only[j-1, i]
                    turbine_l1 += abs(curr_yaw - prev_yaw)
                end
            end
            
            
            # Check Hard Limit for the current turbine
            if turbine_l1 > set_lambda_l1_hard_limit
                return 1e6 + (turbine_l1 - set_lambda_l1_hard_limit)^2 * 1000.0
            end
            
            total_l1_acc += turbine_l1
        end

        # 3. Handle Misalignment Penalty
        if max_violation > 0
            return 1e6 + (max_violation^2 * 1000.0)
        end

        # 4. Final soft penalty
        penalty_term = set_lambda_l1 * total_l1_acc / (wf.nT * (sim.end_time - sim.start_time + 1))


        #otherwise
        try
            state.con.yaw_data = state.yaws_with_time_buffer
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


function construct_yaw_matrix_dynamic_individual_switches!(buffer, x, sim, wf, num_yaw_changes, max_yaw_rate)
    num_steps = size(buffer, 1) # 1201
    duration = sim.end_time - sim.start_time
    num_ts_params = num_yaw_changes - 1
    
    # Calculate the global offset where the yaw values begin
    yaw_start_offset = num_ts_params * wf.nT

    # First Column: Fill with time stamps (sim.start_time : sim.end_time)
    @inbounds for j in 1:num_steps
        buffer[j, 1] = sim.start_time + (j - 1)
    end

    for j in 1:wf.nT
        # Calculate offsets once per turbine
        time_block_start = (j - 1) * num_ts_params
        yaw_block_start = yaw_start_offset + (j - 1) * num_yaw_changes
        
        current_row = 1
        col_idx = j + 1 # +1 because column 1 is time stamps

        for i in 1:num_yaw_changes
            # 1. Determine end_row for this segment
            if i < num_yaw_changes
                # Convert normalized time (0-1) to integer index
                t_val = x[time_block_start + i]
                end_row = round(Int, t_val * duration) + 1
            else
                end_row = num_steps
            end
            
            # Robustness: ensure segments don't overlap backward
            end_row = clamp(end_row, current_row, num_steps)
            
            # 2. Extract yaw value and fill buffer
            val = x[yaw_block_start + i] * 360.0
            
            # Use a simple loop for zero-allocation filling
            @inbounds for row in current_row:end_row
                buffer[row, col_idx] = val
            end
            
            current_row = end_row + 1
            if current_row > num_steps break end
        end
    end

    # 3. Enforce yaw rate limits in-place
    # Note: Ensure apply_yaw_rate_limit! works on a SubArray/View
    yaws_view = @view buffer[:, 2:end]
    apply_yaw_rate_limit!(yaws_view, max_yaw_rate)
    
    return nothing
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


@inline function l1_norm_calc(yaws)
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







#old: --------------------------------------------------------------------------------------




#old implementation for static yaw control
function construct_yaw_matrix(x, sim, wf)   #wf to be used later


    yaws =  ones(sim.end_time - sim.start_time + 1)   * x' * 360.0  #expands into a matrix,
                                                                    # x is in [0,1] range
        
    #ADD CHECK: yaws must not excees X degrees misalignment from wind to prevent crash
    return   [sim.start_time:sim.end_time    yaws]
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



@inline function l1_norm_penalty_slow(yaws)
    # Calculate the L1 norm of the yaw angles, currently farm total
    return sum(abs.(yaws[1:end-1, :] - yaws[2:end, :] )) / wf.nT / (sim.end_time - sim.start_time + 1)  # Normalize by number of turbines and time
end

#axial induction control

function construct_axial_induction_matrix_dynamic!(buffer, x, sim, wf, num_a_changes)
    num_steps = size(buffer, 1)
    duration = sim.end_time - sim.start_time
    
    # 1. Fill the first column with timestamps (in-place)
    @inbounds for j in 1:num_steps
        buffer[j, 1] = sim.start_time + (j - 1)
    end

    # Handle Constant Axial Induction (1 segment)
    if num_a_changes == 1
        for j in 1:wf.nT
            # Assuming x[j] is the induction factor (e.g., 0.33)
            val = x[j] 
            col_idx = j + 1
            @inbounds for i in 1:num_steps
                buffer[i, col_idx] = val
            end
        end
        return nothing
    end

    # 2. Process Dynamic Segments
    current_row = 1
    for i in 1:num_a_changes
        # Determine end_row for this time segment
        if i < num_a_changes
            # x[i] contains the normalized timestamp (0.0 to 1.0)
            end_row = round(Int, x[i] * duration) + 1
        else
            end_row = num_steps
        end
        
        # Safeguard row range
        end_row = clamp(end_row, current_row, num_steps)
        
        if current_row <= end_row
            for j in 1:wf.nT
                # Indexing: (Time params) + (segment_idx * nT) + turbine_idx
                flat_idx = (num_a_changes - 1) + (i - 1) * wf.nT + j
                val = x[flat_idx]
                
                col_idx = j + 1
                @inbounds for row_idx in current_row:end_row
                    buffer[row_idx, col_idx] = val
                end
            end
        end
        current_row = end_row + 1
        if current_row > num_steps break end
    end

    return nothing
end


#main cost function
function test_costfunction_aic(plt, set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, vis::Vis, floridyn::FloriDyn, floris::Floris, objective)
    
    return function cost(x)
        state = tlv[]
        penalty_term = 0.0
    
        #to test AIC: make yaw matrix construction static
        static_yaw_matrix = generate_initial_guess(sim, wind, wf, set_num_yaw_changes)  #generate a static yaw matrix for testing purposes, to isolate the effect of axial induction control on the cost function and gradients 
         
        #yaw matrix construction
        construct_yaw_matrix_dynamic!(state.con.yaw_data, static_yaw_matrix, state.sim, state.wf, set_num_yaw_changes, set_max_yaw_rate)
        
        #get data for penalty function
        sim_times = @view state.con.yaw_data[:, 1]
        yaws_only = @view state.con.yaw_data[:, 2:end] # All turbines   
        
        #get full wind direction matrix
        fill_wind_dir_buffer!(state.wind_dirs_buffer, sim_times, state.wind.dir)


        #generate axial induction matrix
        #@infiltrate
        construct_axial_induction_matrix_dynamic!(state.con.induction_data, x, state.sim, state.wf, set_num_a_changes)
        

        #handle soft and hard constraints
        penalty_term = calculate_penalties(yaws_only, state.wind_dirs_buffer)
        
        #detect if penalty term exceeds hard limit, and if so return it directly to avoid unnecessary simulations
        if penalty_term >= 1e6
            return penalty_term
        end
        
        #otherwise
        try
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


function generate_initial_guess_aic(sim, wind, wf, n_segments)
    # Handle the constant induction case (1 segment)
    if n_segments == 1
        return fill(set_initial_a_value, wf.nT)
    end

    # 1. Internal transition times (normalized 0 to 1)
    # This creates [0.0, 0.25, 0.5, 0.75, 1.0] for n_segments = 4
    equal_time_spacing = collect(0:1/n_segments:1)
    
    # Extract only the internal knots (e.g., [0.25, 0.5, 0.75])
    time_guesses = equal_time_spacing[2:end-1]
    
    # 2. Induction factor guesses
    # We need a value for every turbine (nT) in every time segment (n_segments)
    num_vals = n_segments * wf.nT
    a_guesses = fill(Float64(set_initial_a_value), num_vals)
    
    # Concatenate: [timestamps..., induction_values...]
    return vcat(time_guesses, a_guesses)
end




#testing area inside testing area----------------------------------------------------------------------------
