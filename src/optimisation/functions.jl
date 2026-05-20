using Base.Threads



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
    nrows, ncols = size(yaws)
    @inbounds for i in 2:nrows
        @inbounds for j in 1:ncols
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


#fill wind direction buffer
function fill_wind_dir_buffer!(buffer, sim_times, wind_matrix)
    # dest_buffer: N elements
    # sim_times: N elements
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
        
        state = tlv[][]
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
            return opt_set.set_objective(state.power_vector, state.wf.nT, state.sim.n_sim_steps, opt_set.set_num_timesteps_to_skip) + penalty_term  #in kW per turbine

        catch e
            # 1. If the user hits Ctrl+C, let it happen!
            if e isa InterruptException 
                rethrow(e)
            end
            if e isa ArgumentError #memory errors for Evolutionary.jl fall under this
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
@inline function totalEnergyObjective(powervector, nT, nsteps, num_timesteps_to_skip=0)

    #return -sum(powervector) / (nT * nsteps) * 1000.0   #in kW per turbine
    
    #num_timesteps_to_skip is to skip the first 'uninitialised' timesteps
    #value is given by largest_length_scale / u / dt
    #@infiltrate
        
        return -sum(@view powervector[(nT*num_timesteps_to_skip+1):end] ) / (nT * (nsteps - num_timesteps_to_skip)) * 1000.0
        
        
end


function totalEnergyObjectivePlusCurvatureRegulariser(
    powervector, nT, nsteps, num_timesteps_to_skip=125, lambda=3e5
)
    active_start = nT * num_timesteps_to_skip + 1
    nActive_steps = nsteps - num_timesteps_to_skip

    # Primary objective: maximise mean power (same as before)
    energy_sum = zero(eltype(powervector))
    for i in active_start:length(powervector)
        energy_sum += powervector[i]
    end
    energy_term = -energy_sum / (nT * nActive_steps) * 1000.0

    # Curvature regulariser: discrete second differences on aggregate power per timestep
    curvature = zero(eltype(powervector))
    @inbounds for t in 2:(nActive_steps - 1)
        p_prev = zero(eltype(powervector))
        p_curr = zero(eltype(powervector))
        p_next = zero(eltype(powervector))
        @inbounds for j in 1:nT
            p_prev += powervector[active_start + nT*(t-2) + j - 1]
            p_curr += powervector[active_start + nT*(t-1) + j - 1]
            p_next += powervector[active_start + nT*(t-0) + j - 1]
        end
        curvature += (p_next - 2*p_curr + p_prev)^2
    end
    curvature_term = lambda * curvature / (nActive_steps - 2)

    return energy_term + curvature_term
end




#wind direction interpolation helper function
#in-loop replaced by fill_wind_dir_buffer! for performance
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






