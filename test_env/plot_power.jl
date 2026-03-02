#plotting the power:



using Plots




#set_desired_power_curve
function get_power_curve_from_measurements(md,sim)


    power_per_step = zeros(sim.n_sim_steps)
    for k in 1:sim.n_sim_steps
    acc = 0.0
    offset = (k - 1) * wf.nT
    @inbounds for j in 1:wf.nT
        acc += md.PowerGen[offset + j]
            end
    power_per_step[k] = acc
    end
    
    return power_per_step# in MW
end


optimised_power_curve = get_power_curve_from_measurements(md, sim)


#plot both set_desired_power_curve and optimised_power_curve



using Plots
plot(1:sim.n_sim_steps, set_desired_power_curve, label="desired power curve", xlabel="time step", ylabel="power (kW)")
plot!(1:sim.n_sim_steps, optimised_power_curve, label="optimised power curve")






