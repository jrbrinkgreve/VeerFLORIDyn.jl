using Plots
plot = Plots.plot
savefig = Plots.savefig
#plotting the regularisation results

function regularisation_results()
    #note: 100 is surrogate for 0 here
        lambda_vals = [500, 1000, 2000, 3000, 5000, 7000, 10000, 20000, 30000, 50000]
        lambda_labels = ['0', "1e3", "2e3", "3e3", "5e3", "7e3", "1e4", "2e4", "3e4", "5e4"]

        # Increase over baseline (%)
        P_pct = [1.72, 1.70, 1.71, 1.65, 1.5, 1.26, 0.85, 0.0, -0.75, -1.74]

        # Optimized average L1 yaw change norm
        norm1 = [51.91, 49.95, 48.31, 44.69, 39.25, 32.78, 24.77, 14.72, 9.83, 5.14]

    p1 = plot(lambda_vals, norm1,
        xscale=:log10,
        ylabel="‖Δγ_opt‖₁\n(°)",
        ylims=(0.0, 60.0),
        yticks=0.0:5.0:60.0,
        label="‖Δγ_opt‖₁",
        color=:coral,
        marker=:diamond,
        linewidth=2,
        title="Regularised Optimisation Results",
        xticks=(lambda_vals, lambda_labels),
        xlabel="λ",
        legendfontsize=8,
        legend=:bottomleft,
        grid=true,
        gridalpha=0.4)
        

    plot!(twinx(), lambda_vals, P_pct,
        xscale=:log10,
        ylabel="P%\n(%)",
        ylims=(-3.0, 3.0),
        yticks=-3.0:0.5:3.0,
        label="P%",
        color=:mediumseagreen,
        marker=:circle,
        linewidth=2,
        legend=:topright,
        grid = false)

    plot(p1,
        size=(700, 350),
        left_margin=5Plots.mm,
        right_margin=15Plots.mm,
        ygrid=true,
        xgrid = true,
        gridalpha=0.4,
        bottom_margin=5Plots.mm)

    savefig("output/regularisation_results.pdf")
end





#------------------------------------------------------------------------------------------------------------------
#plotting the convexity analysis grid search results

using JLD2


function convexity_analysis_results()
    # Load the matrices


    @load "src/optimisation/optim_data/matrix_noveer.jld2" matrix_noveer
    @load "src/optimisation/optim_data/matrix_0.01veer.jld2" matrix_0_01veer
    @load "src/optimisation/optim_data/matrix_0.05veer.jld2" matrix_0_05veer
    @load "src/optimisation/optim_data/matrix_0.1veer.jld2" matrix_0_1veer
    @load "src/optimisation/optim_data/matrix_0.15veer.jld2" matrix_0_15veer
    @load "src/optimisation/optim_data/matrix_0.2veer.jld2" matrix_0_2veer


    N = size(matrix_noveer, 1)  # Assuming all matrices have the same dimensions
    yaw_lims = 80.0 #degrees, the range of yaw misalignments considered in the grid search

    # Define the corresponding wind veer gradients for each matrix
    veer_gradients = [0.0, 0.01, 0.05, 0.1, 0.15, 0.2]

   
end




function l0_norm_results()
    norm0_vals = [2, 3, 4, 5, 6]
    P_opt  = [3910.74, 3938.71, 3983.64, 3989.31, 3996.65]
    P_pct  = [-0.16, 0.57, 1.72, 1.86, 2.05]
    norm1  = [20.45, 28.33, 51.91, 55.87, 61.11]

    p1 = plot(norm0_vals, norm1,
        ylabel="‖Δγ_opt‖₁\n(°)",
        ylims=(10.0, 70.0),
        yticks=10.0:10.0:70.0,
        label="‖Δγ_opt‖₁",
        color=:coral,
        marker=:diamond,
        linewidth=2,
        title="L0 Norm Constraint Results",
        xticks=(norm0_vals, string.(norm0_vals)),
        xlabel="‖Δγ‖₀  (No. turbine actuations)",
        legendfontsize=8,
        legend=:topleft,
        grid=true,
        gridalpha=0.4)

        
    plot!(twinx(), norm0_vals, P_pct,
        ylabel="P%\n(%)",
        ylims=(-0.5, 2.5),
        yticks=-0.5:0.5:2.5,
        label="P%",
        color=:mediumseagreen,
        marker=:circle,
        linewidth=2,
        legend=:bottomright,
        grid=false)

    plot(p1,
        size=(700, 350),
        left_margin=5Plots.mm,
        right_margin=15Plots.mm,
        ygrid=true,
        xgrid=true,
        gridalpha=0.4,
        bottom_margin=5Plots.mm)

    savefig("output/l0_norm_results.pdf")
end


function veer_results()
    alpha_vals = [0.0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2]

    P_base = [3916.40, 3914.99, 3904.37, 3884.84, 3859.09, 3814.67, 3740.25, 3678.85]
    P_opt  = [3983.64, 3981.16, 3966.13, 3944.19, 3914.75, 3866.63, 3790.86, 3728.18]
    P_pct  = [1.72, 1.69, 1.58, 1.53, 1.44, 1.36, 1.35, 1.34]
    norm1  = [51.91, 51.46, 51.55, 50.27, 49.21, 46.51, 41.85, 39.49]

    alpha_labels = string.(alpha_vals)

   p1 = plot(alpha_vals, P_base,
        ylabel="Power\n(kW)",
        ylims=(3650.0, 4050.0),
        yticks=3650.0:50.0:4050.0,
        label="P_base",
        color=:steelblue,
        marker=:circle,
        linewidth=2,
        title="Optimal Yaw Steering Gain vs. Veer Gradient Strength",
        xticks=(alpha_vals, alpha_labels),
        xlabel="α (deg/m)",
        legendfontsize=8,
        legend=:bottomleft,
        xrot=45,
        grid=true,
        gridalpha=0.4)

    plot!(alpha_vals, P_opt,
        label="P_opt",
        color=:mediumpurple,
        marker=:circle,
        linewidth=2)

    plot!(twinx(), alpha_vals, P_pct,
        ylabel="P%\n(%)",
        ylims=(1.0, 2.0),
        yticks=1.0:0.25:2.0,
        label="P%",
        color=:mediumseagreen,
        marker=:diamond,
        linewidth=2,
        legendfontsize=8,
        legend=:topright,
        grid=false)

    plot(p1,
        size=(700, 350),
        left_margin=5Plots.mm,
        right_margin=15Plots.mm,
        ygrid=true,
        xgrid=true,
        gridalpha=0.4,
        bottom_margin=5Plots.mm)

    savefig("output/veer_results.pdf")
end





function yaw_limit_results()
    # Yaw limit values (using 999 as surrogate for Inf)
    yaw_limits     = [10, 15, 20, 25, 30]
    yaw_labels     = ["10°", "15°", "20°", "25°", "Inf"]

    # Increase over baseline (%)
    P_pct  = [0.55, 1.11, 1.49, 1.72, 1.77]

    # Optimized average L1 yaw change norm (deg)
    norm1  = [31.27, 36.92, 44.09, 51.44, 54.79]

    p1 = Plots.plot(yaw_limits, norm1,
        ylabel="‖Δγ_opt‖₁\n(°)",
        ylims=(15.0, 60.0),
        yticks=15.0:5:60.0,
        label="‖Δγ_opt‖₁",
        color=:coral,
        marker=:diamond,
        linewidth=2,
        title="Yaw Misalignment Limit Optimisation Results",
        xticks=(yaw_limits, yaw_labels),
        xlabel="Yaw Misalignment Limit (°)",
        legendfontsize=8,
        legend=:topleft,
        grid=true,
        gridalpha=0.4)
    Plots.plot!(twinx(), yaw_limits, P_pct,
        ylabel="P%\n(%)",
        ylims=(0.0, 2.25),
        yticks=0.0:0.25:2.25,
        label="P%",
        color=:mediumseagreen,
        marker=:circle,
        linewidth=2,
        legend=:bottomright,
        grid=false)
    Plots.plot(p1,
        size=(700, 350),
        left_margin=5Plots.mm,
        right_margin=15Plots.mm,
        ygrid=true,
        xgrid=true,
        gridalpha=0.4,
        bottom_margin=5Plots.mm)
    Plots.savefig("output/yaw_limits_results.pdf")
end












#regularisation_results()
#convexity_analysis_results()
l0_norm_results()
#veer_results()
#runtime_results_threads()
yaw_limit_results()




nothing