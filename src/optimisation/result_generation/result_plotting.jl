using Plots
#plotting the regularisation results
function regularisation_results()
    lambda_vals   = [1e3, 5e3, 1e4, 2e4, 5e4]
    lambda_labels = ["10³", "5×10³", "10⁴", "2×10⁴", "5×10⁴"]

    P_pct  = [2.51,  2.44,  1.71,  0.25,  -1.35]
    norm1  = [51.01, 46.93, 32.60, 15.86,   6.58]

    p1 = plot(lambda_vals, norm1,
        xscale=:log10,
        ylabel="‖Δγ_opt‖₁\n(°)",
        ylims=(0.0, 65.0),
        yticks=0.0:10.0:65.0,
        label="‖Δγ_opt‖₁",
        color=:coral,
        marker=:diamond,
        linewidth=2,
        title="Regularised Optimisation Results",
        xticks=(lambda_vals, lambda_labels),
        xlabel="λ",
        legendfontsize=8,
        legend=:topright,
        gridalpha=0.4)

    plot!(twinx(), lambda_vals, P_pct,
        xscale=:log10,
        ylabel="P%\n(%)",
        ylims=(-2.0, 3.0),
        yticks=-2.0:0.5:3.0,
        label="P%",
        color=:mediumseagreen,
        marker=:circle,
        linewidth=2,
        legend=:bottomleft,
        grid=false)

    plot(p1,
        size=(700, 350),
        left_margin=5Plots.mm,
        right_margin=15Plots.mm,
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

        P_opt  = [3859.30, 3883.99, 3951.13, 3962.80, 3964.25]
        P_pct  = [0.11,    0.75,    2.49,    2.80,    2.83]
        norm1  = [22.24,   28.32,   54.18,   58.65,   60.79]

        #=
        p1 = plot(norm0_vals, P_opt,
            ylabel="P_opt (kW)",
            label="P_opt",
            color=:steelblue,
            marker=:circle,
            linewidth=2,
            title="Optimal Power",
            xticks=(norm0_vals, string.(norm0_vals)),
            legend=:topleft)
        =#
        p2 = plot(norm0_vals, P_pct,
            ylabel="P% (%)",
            label="P% with ‖Δγ_base‖₀ = 4",


            color=:mediumseagreen,
            marker=:circle,
            linewidth=2,
            title="Relative Power increase & Control Effort Magnitude",
            xticks=(norm0_vals, string.(norm0_vals)),
            xlabel="‖Δγ‖₀",
            legendfontsize=8,
            legend=:topleft,
    
            yticks=0:0.5:3.0,
    
            gridalpha=0.4,
            ylims=(0, 3.0))


        plot!(twinx(), norm0_vals, norm1,
            ylabel="‖Δγ_opt‖₁ (°)",
            label="‖Δγ_opt‖₁",
            color=:coral,
            marker=:diamond,
            linewidth=2,
            ylims=(0, 70),
            grid=false,
            legend=:bottomright)

        plot(p2,
            size=(600, 300),
        )


    savefig("output/l0_norm_results.pdf")

end



function veer_results()
    alpha_vals = [0, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2]

    P_base = [3855.04, 3854.72, 3853.76, 3825.84, 3764.04, 3702.19, 3652.95]
    P_opt  = [3951.13, 3949.78, 3947.03, 3912.35, 3841.17, 3768.88, 3717.84]
    P_pct  = [2.49,    2.47,    2.42,    2.26,    2.05,    1.80,    1.78]
    norm1  = [54.18,   52.75,   53.64,   52.36,   49.25,   45.58,   42.51]

    alpha_labels = string.(alpha_vals)

    p1 = plot(alpha_vals, P_base,
        ylabel="Power (kW)",
        ylims=(3000.0, 4000.0),
        label="P_base",
        color=:steelblue,
        marker=:circle,
        linewidth=2,
        title="Optimal yaw steering gain vs. Veer gradient strength",
        xticks=(alpha_vals, alpha_labels),
        xlabel="α [deg/m]",
        legendfontsize=10,
        legend=:bottomleft,
        xrot=45,
        gridalpha=0.4)

    plot!(alpha_vals, P_opt,
        label="P_opt",
        color=:mediumpurple,
        marker=:circle,
        linewidth=2)

    plot!(twinx(), alpha_vals, P_pct,
        ylabel="P%\n(%)",
        ylims=(1.0, 3.0),
        yticks=1.0:0.5:3.0,
        label="P%",
        color=:mediumseagreen,
        marker=:diamond,
        linewidth=2,
        legendfontsize=10,
        legend=:bottomright,
        grid=false)

    plot(p1,
        size=(800, 400),
        left_margin=5Plots.mm,
        right_margin=15Plots.mm,
        bottom_margin=5Plots.mm)

    savefig("output/veer_results.pdf")
end


regularisation_results()
#convexity_analysis_results()
#l0_norm_results()
veer_results()



nothing