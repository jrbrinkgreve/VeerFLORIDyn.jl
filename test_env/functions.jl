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


#cost function with parallel
function parallel_costfunction(plt, set, wf::WindFarm, wind::Wind, sim, con, vis, floridyn, floris)
    tlv = TaskLocalValue{NamedTuple}() do
        (
            plt=deepcopy(plt), set=deepcopy(set), wf=deepcopy(wf),
            wind=deepcopy(wind), sim=deepcopy(sim), floris=deepcopy(floris),
            floridyn=deepcopy(floridyn), vis=deepcopy(vis), con=deepcopy(con)
        )
    end
    
    return function cost(x)
        state = tlv[]
        
        # 1. Construct the yaw matrix [time y1 y2 ...]
        yaws_with_time = construct_yaw_matrix_dynamic(x, state.sim, state.wf, set_num_yaw_changes, set_max_yaw_rate)
        
        # 2. Extract Data for Validation
        sim_times = yaws_with_time[:, 1]
        yaws_only = yaws_with_time[:, 2:end] # All turbines
        wind_table = state.wind.dir           # Your 3x2 matrix
        
        # 3. Calculate Relative Yaw Violation
        # We check every time step 't' against the interpolated wind at 't'
        max_violation = 0.0
        @inbounds for (i, t) in enumerate(sim_times)
            current_wind_dir = get_wind_at_t(t, wind_table)
          
            # Check all turbines at this specific second
            @inbounds for turbine_idx in axes(yaws_only, 2)
                rel_yaw = abs(yaws_only[i, turbine_idx] - current_wind_dir)
                # Keep track of the worst offender relative to the 90 deg crash limit
                # We use 85 deg as a buffer
                max_violation = max(max_violation, rel_yaw - 85.0)

            end
        end

        # 4. Penalty Branch
        if max_violation > 0
            # Quadratic penalty creates a smooth "slope" leading back to 0 violation
            return 1e6 + (max_violation^2 * 1000.0)
        end

        # 5. Feasible Branch
        try
            state.con.yaw_data = yaws_with_time
            wf, md, mi = run_floridyn(
                state.plt, state.set, state.wf, state.wind, 
                state.sim, state.con, state.vis, state.floridyn, state.floris
            )
            
            
            #can add more complex term here
            #is now set to total energy
            return -sum(md.PowerGen) / (wf.nT * sim.n_sim_steps) * 1000.0  #in kW per turbine




        catch e
            # 1. If the user hits Ctrl+C, let it happen!
            if e isa InterruptException
                rethrow(e)
            end
            return 2e6 # Fallback for unexpected failures
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