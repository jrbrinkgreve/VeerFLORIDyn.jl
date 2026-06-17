using JLD2
using FLORIDyn

include("../src/optimisation/optimisationstructs.jl")
include("../src/optimisation/functions.jl")
include("../src/optimisation/MPC_functions.jl")

artifact_file = joinpath("output", "mpc_testbench", "mpc_artifacts.jld2")
if !isfile(artifact_file)
    error("Missing MPC artifact file: $(artifact_file). Run examples/mpc_yaw_testbench.jl first.")
end

data = JLD2.load(artifact_file, "data")
plot_mpc_visualisation(
    data;
    save_path = joinpath("output", "mpc_testbench", "mpc_yaw_control.pdf"),
    show_plot = true,
    title_prefix = "MPC testbench",
)

println("Saved MPC visualisation to output/mpc_testbench/mpc_yaw_control.pdf")
