#goal: validate results with a full wind farm simulation
using Statistics


include("../../../examples/veer_main_mini.jl")

num_timesteps_to_skip = 125
powervector = md.PowerGen[1+(num_timesteps_to_skip+1)*wf.nT:end] * 1000.0

# Store power time series per turbine
n_steps = length(powervector) ÷ wf.nT
turbine_power_matrix = reshape(powervector, wf.nT, n_steps)'  # shape: (n_steps × nT)



# Compute and print average power per turbine
for t in 1:wf.nT
    avg = mean(turbine_power_matrix[:, t])
    println("Turbine $t average power: $avg kW")
end