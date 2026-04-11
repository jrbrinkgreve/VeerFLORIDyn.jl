
using BenchmarkTools


#yaw and veer power production
eps = 1e-3
yaw = 10.0 # 20.0 #deg
veer_gradient = 0.05 #0.1
z_hub = 119.0
R = 89.2
p_p = 2.2


dz = 0.01



z = (z_hub-R+eps):dz:(z_hub+R-eps)
h_sq = R^2 .- (z .- z_hub).^2
eff_yaw = deg2rad.(yaw .+ veer_gradient .* (z .- z_hub)    )



eff_area = 2.0 * sum(    sqrt.(h_sq) .*  cos.(eff_yaw).^p_p .* dz    )
println("Numerical: ")
println(eff_area)

println()
#analytical approx
eff_area_analytical = pi * R^2 * cosd(yaw)^p_p * (1.0 - (p_p * deg2rad(veer_gradient)^2 * R^2 ) / 8.0 *  (1.0 -(p_p-1.0) * tand(yaw)^2 ) ) 

println("Analytical approx: ")
println(eff_area_analytical)




#power production:

#P = 0.5 * rho * A_eff * v^3 * Cp

Cp = 0.45
rho = 1.225
v = 8.0
P = 0.5 * rho * pi*R^2 * v^3 * Cp * cosd(yaw)^p_p
println()
println("Power production with yaw and veer: ")
println(P, " W")
println()
println("with correction:")
P = 0.5 * eff_area_analytical * v^3 * Cp * rho
println(P, " W")
