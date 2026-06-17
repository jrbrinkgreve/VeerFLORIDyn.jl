#a file with all MPC functions








function load_all_data()
    #get some random wind direction, speed, turbulence and veer data for the MPC to work with
end

function write_horizon_limited_csvs(horizon, time, data, path_to_CSVs)
    #update the CSVs in data folder for the optimiser to read, e.g. for floridyn to use
end



function parallel_yaw_optimisation_MPC(init)
    #this is the big function, run the optimisation loop but with a different 
end

function display_MPC_control()
    
end










#main cost function
function parallel_costfunction_MPC(plt, set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, vis::Vis, floridyn::FloriDyn, floris::Floris, opt_set::OptimisationSettings)
    
    return function cost(x)
        
        state = tlv[][]
        penalty_term = 0.0
    
        #x is size (num_yaw_changes - 1 + (num_yaw_changes - 1) * wf.nT)
        


        #yaw matrix construction
        construct_yaw_matrix_MPC_dynamic!(state.con.yaw_data, x, prev_x, t,  state.sim, state.wf, opt_set)
        
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









function construct_yaw_matrix_MPC_dynamic!(buffer, x, prev_x, t,  sim, wf, opt_set)
    num_yaw_changes = opt_set.set_num_yaw_changes
    max_yaw_rate = opt_set.set_max_yaw_rate
    num_steps = size(buffer, 1)
    duration = sim.end_time - sim.start_time

    # 1. Fill the first column with timestamps (in-place)
    @inbounds for j in 1:num_steps
        buffer[j, 1] = sim.start_time + (j - 1)
    end

    # 2. Fill the fixed first segment from prev_x, up to row corresponding to t
    t_row = round(Int, t - sim.start_time) + 1
    t_row = clamp(t_row, 1, num_steps)

    for j in 1:wf.nT
        # Extract yaw values from prev_x using the same indexing as the original constructor
        flat_idx = (num_yaw_changes - 1) + j   # i=1 segment
        val = prev_x[flat_idx] * 360.0
        col_idx = j + 1
        @inbounds for row_idx in 1:t_row
            buffer[row_idx, col_idx] = val
        end
    end

    # 3. Process remaining segments from x
    # x layout: [t_valid, t_3, ..., t_{m-1}, γ_21..γ_2n, γ_31..γ_3n, ...]
    # Number of free time params = num_yaw_changes - 2 (t_valid + interior timestamps)
    num_free_time_params = num_yaw_changes - 1  # t_valid, t_3, ..., t_{m-1}
    current_row = t_row + 1

    for i in 2:num_yaw_changes
        # Determine end_row for this segment
        if i < num_yaw_changes
            # x[i-1] is the normalised timestamp for this segment
            end_row = round(Int, x[i - 1] * duration) + 1
        else
            end_row = num_steps
        end

        end_row = clamp(end_row, current_row, num_steps)

        if current_row <= end_row
            for j in 1:wf.nT
                # Yaw values start after the free time params in x
                flat_idx = num_free_time_params + (i - 2) * wf.nT + j
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

    # 4. Enforce yaw rate limits in-place on the yaw columns
    apply_yaw_rate_limit!(@view(buffer[:, 2:end]), max_yaw_rate)
    return nothing
end
