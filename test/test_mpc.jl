using Test
using FLORIDyn

include("../src/optimisation/optimisationstructs.jl")
include("../src/optimisation/functions.jl")
include("../src/optimisation/MPC_functions.jl")

@testset "MPC helpers" begin
    @testset "time window limiting" begin
        data = [0.0 270.0 8.0; 10.0 280.0 9.0]
        limited = _limit_time_series(data, 2.0, 8.0)

        @test size(limited, 2) == 3
        @test limited[1, 1] ≈ 2.0
        @test limited[end, 1] ≈ 8.0
        @test limited[1, 2] ≈ 272.0
        @test limited[end, 2] ≈ 278.0
    end

    @testset "visualisation output" begin
        time = collect(0.0:1.0:4.0)
        yaw_matrix = hcat(fill(270.0, length(time)), fill(275.0, length(time)))
        wind_dirs = fill(272.0, length(time))
        output_dir = mktempdir()
        output_file = joinpath(output_dir, "mpc_yaw_control.pdf")

        p = plot_mpc_visualisation(time, yaw_matrix, wind_dirs; save_path = output_file, show_plot = false, title_prefix = "MPC smoke test")

        @test isfile(output_file)
        @test p !== nothing
    end
end
