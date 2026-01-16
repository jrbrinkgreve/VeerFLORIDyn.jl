

#yaw and veer power production
eps = 1e-3
yaw = 0.0 # 20.0 #deg
veer_gradient = 0.05 #0.1
z_hub = 119.0
R = 89.2
dz = 0.000001



z = (z_hub-R+eps):dz:(z_hub+R-eps)
h_sq = R^2 .- (z .- z_hub).^2
eff_yaw = deg2rad.(yaw .+ veer_gradient .* (z .- z_hub)    )
p_p = 2.2


eff_area = 2.0 * sum(    sqrt.(h_sq) .*  cos.(eff_yaw).^p_p .* dz    )







        