
using BenchmarkTools
using Printf

#yaw and veer power production


#=

eps = 1e-3
yaw = 30.0 # 20.0 #deg
veer_gradient = 0.2 #0.1
z_hub = 119.0
R = 89.2
p_p = 2.2
dz = 1e-5

z = (z_hub-R+eps):dz:(z_hub+R-eps)
h_sq = R^2 .- (z .- z_hub).^2
eff_yaw = deg2rad.(yaw .+ veer_gradient .* (z .- z_hub)    )



eff_numerical = 2.0 * sum(    sqrt.(h_sq) .*  cos.(eff_yaw).^p_p .* dz    ) / (pi * R^2)
println("Numerical: ")
println(eff_numerical)

println()
#analytical approx
eff_analytical =  cosd(yaw)^p_p * (1.0 - (p_p * deg2rad(veer_gradient)^2 * R^2 ) / 8.0 *  (1.0 -(p_p-1.0) * tand(yaw)^2 ) ) 

println("Analytical approx: ")
println(eff_analytical)


error = abs(eff_numerical - eff_analytical) / abs(eff_numerical)
println()
println("Relative error: ", error * 100, " %")  


=#


# Yaw and veer power production - parameter sweep
# Computes numerical vs analytical efficiency and relative error
# for veer in 0:0.05:0.2 and yaw in 0:5:30

eps     = 1e-3
z_hub   = 119.0
R       = 89.2
p_p     = 2.2
dz      = 1e-5

veer_range = 0.0:0.05:0.2
yaw_range  = 0.0:5.0:30.0

# Header
println("Veer  | Yaw   | Numerical  | Analytical | Rel. Error (%)")
println(repeat("-", 62))

for veer_gradient in veer_range
    for yaw in yaw_range
        z = (z_hub - R + eps):dz:(z_hub + R - eps)
        h_sq = R^2 .- (z .- z_hub).^2

        eff_yaw       = deg2rad.(yaw .+ veer_gradient .* (z .- z_hub))
        eff_numerical = 2.0 * sum(sqrt.(h_sq) .* cos.(eff_yaw).^p_p .* dz) / (pi * R^2)

        eff_analytical = cosd(yaw)^p_p * (
            1.0 - (p_p * deg2rad(veer_gradient)^2 * R^2) / 8.0 *
            (1.0 - (p_p - 1.0) * tand(yaw)^2)
        )

        rel_error = abs(eff_numerical - eff_analytical) / abs(eff_numerical) * 100

        @printf("%-6.2f| %-6.1f| %-11.6f| %-11.6f| %.4f\n",
                veer_gradient, yaw, eff_numerical, eff_analytical, rel_error)
    end
    println(repeat("-", 62))
end


#power production:

#P = 0.5 * rho * A_eff * v^3 * Cp
#=

Cp = 0.45
rho = 1.225
v = 8.0
P = 0.5 * rho * pi*R^2 * v^3 * Cp * cosd(yaw)^p_p
println()
println("Power production with yaw and veer: ")
println(P, " W")
println()
println("with correction:")
P = 0.5 * eff_analytical * v^3 * Cp * rho
println(P, " W")

=#