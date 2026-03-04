@kwdef struct OptimisationSettings
    set_num_yaw_changes::Int  #N
    set_num_optimiser_runs::Int#number of automatic restarts for CMA-ES
    set_max_yaw_rate::Float64 #deg/s
    set_max_yaw_misalignment::Float64 
    set_lambda_l1::Float64 
    set_lambda_l1_hard_limit::Float64
    set_objective::Function 
end



