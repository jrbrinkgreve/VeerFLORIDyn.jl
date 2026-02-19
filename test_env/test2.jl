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
