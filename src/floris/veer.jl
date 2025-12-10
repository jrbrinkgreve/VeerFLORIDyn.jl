#this file includes veer calculation functions


#=
List of functions in this file:

- calcCT function already  in gaussian.jl, not defined here due to double method defining errors
- getVars! function for getting veer related variables: added dispatch method with veer_enabled == 1
- centerline! function: TBD what to do with this
- mutable struct States?
- constructor for States?
- init_states: probably just copy this from gaussian.jl
- getPower: function to get power, with veer 
- getUadv: function to get Uadv, with veer

=#




"""

                             BUFFERS / INPLACE           |||        INPUT VARS
        getVars!(sig_y, sig_z, x_0, delta, pc_y, pc_z, rps, c_t, yaw, ti, ti0, floris::Floris, d_rotor, veer_enabled)

Compute Gaussian wake widths, deflection, potential-core radii, and onset distance at observation points, in-place.

# Output Arguments
- `sig_y::AbstractVector{<:Real}` (length n): Lateral Gaussian width σ_y at each point [m]
- `sig_z::AbstractVector{<:Real}` (length n): Vertical Gaussian width σ_z at each point [m]
- `x_0::AbstractVector{<:Real}` (length n): Onset distance of the far-wake x₀ [m]
- `delta::AbstractMatrix{<:Real}` (length n×2): Deflection components `[Δy, Δz]` [m]
- `pc_y::AbstractVector{<:Real}` (length n): Potential-core radius in y at each point [m]
- `pc_z::AbstractVector{<:Real}` (length n): Potential-core radius in z at each point [m]





# Input Arguments
- `rps::AbstractMatrix` (n×3): Observation points in wake-aligned frame; columns are `[x_downstream, y_cross, z_cross]` [m]
- `c_t::Union{Number,AbstractVector}`: Thrust coefficient Ct (scalar or length n) [-]
- `yaw::Union{Number,AbstractVector}`: Yaw misalignment (scalar or length n) [rad]
- `ti::Union{Number,AbstractVector}`: Local turbulence intensity TI at turbine (scalar or length n) [-]
- `ti0::Union{Number,AbstractVector}`: Ambient turbulence intensity TI₀ (scalar or length n) [-]
- `floris::Floris`: FLORIS Gaussian model parameters; see [`Floris`](@ref)
- `d_rotor::Real`: Rotor diameter D [m]


"""

function getVars!(sig_y, sig_z, x0, delta, pc_y, pc_z, RPs, Ct, yaw, TI, TI0, floris, D, veer_enabled::Int)
    #note: wake aligned frame here means at hub height, and per definition veer angle at yaw angle is 0.0 deg
    

    


    return nothing
end




#=
use 2023 wake combination paper by Li et al. for wake combination using sum of squares for velocity deficits



=#
