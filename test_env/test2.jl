

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
    for i in 1:(num_yaw_changes - 1)
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
            for j in 1:wf.nT
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
