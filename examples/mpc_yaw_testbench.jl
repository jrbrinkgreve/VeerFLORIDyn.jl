using FLORIDyn

include("../src/optimisation/optimisationstructs.jl")
include("../src/optimisation/functions.jl")
include("../src/optimisation/MPC_functions.jl")

ctx, data = run_mpc_testbench(
    vis_file = "data/vis_default.yaml",
    settings_file = "data/MPC_formulation.yaml",
    output_dir = "output/mpc_testbench",
    horizon = 120,
    receding_step = 60,
    num_yaw_changes = 2,
    num_optimiser_runs = 1,
    iterations = 3,
    num_timesteps_to_skip = 0,
    sigma0 = 0.03,
    show_plot = false,
)

println("MPC testbench completed")
println("Objective: $(data.objective)")
println("Artifacts written to output/mpc_testbench")
