# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    FLORISBuffers

Pre-allocated buffers for the runFLORIS! computation to minimize allocations.

This struct also persists result arrays so callers can read outputs without
allocations. After calling `runFLORIS!`, the following fields contain results:

# Output Fields
- `T_red_arr::Vector{Float64}`: Per-turbine velocity reduction factors. For a
    single-turbine run, length is 1 and `T_red_arr[1]` is the scalar reduction.
- `T_aTI_arr::Vector{Float64}`: Added turbulence intensity from upstream wakes.
    For N turbines, length is `max(N-1, 0)`. Empty for single-turbine runs.
- `T_Ueff::Vector{Float64}`: Effective wind speed at the last turbine as a
    length-1 vector (multi-turbine case). Empty for single-turbine runs.
- `T_weight::Vector{Float64}`: Gaussian weight factors used for wake overlap.
    For N turbines, length is `max(N-1, 0)`. Empty for single-turbine runs.
"""
mutable struct FLORISBuffers
    tmp_RPs::Matrix{Float64}
    rotor_pts::Matrix{Float64}
    # Preallocated arrays for getVars! outputs
    sig_y::Vector{Float64}
    sig_z::Vector{Float64}
    x_0::Vector{Float64}
    delta::Matrix{Float64}   # n×2
    pc_y::Vector{Float64}
    pc_z::Vector{Float64}
    cw_y::Vector{Float64}
    cw_z::Vector{Float64}
    phi_cw::Vector{Float64}
    r_cw::Vector{Float64}
    core::Vector{Bool}
    nw::Vector{Bool}
    fw::Vector{Bool}
    tmp_RPs_r::Vector{Float64}
    gaussAbs::Vector{Float64}
    gaussWght::Vector{Float64}
    exp_y::Vector{Float64}
    exp_z::Vector{Float64}
    not_core::Vector{Bool}



    
    nRP::Int
    rps_coords::Matrix{Float64}  # 3 x nRP
    alpha::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    a_star::Vector{Float64}
    xi_0_hat::Vector{Float64}
    rotmtx::Matrix{Float64}
    coords_veered::Matrix{Float64}
    shear_modifier::Vector{Float64}
    u_in_z::Vector{Float64}
    t_hat::Vector{Float64}
    sgn_t_hat::Vector{Float64}
    abs_t_hat::Vector{Float64}
    y_hat_c::Vector{Float64}
    y_c::Vector{Float64}
    theta::Vector{Float64}
    xi_0::Vector{Float64}
    xi_hat::Vector{Float64}
    chi::Vector{Float64}
    a::Vector{Float64}
    c1::Vector{Float64}
    c2::Vector{Float64}
    c3::Vector{Float64}
    c4::Vector{Float64}
    c5::Vector{Float64}
    c6::Vector{Float64}
    c7::Vector{Float64}
    xi::Vector{Float64}
    sigma::Vector{Float64}
    sigma_hat_squared::Vector{Float64}
    c::Vector{Float64}
    du::Vector{Float64}
    u::Vector{Float64}
    tmp::Vector{Float64}

    # Result arrays (persisted in buffers to avoid fresh allocations)
    T_red_arr::Vector{Float64}
    T_aTI_arr::Vector{Float64}
    T_Ueff::Vector{Float64}    # length 1 when set
    T_weight::Vector{Float64}
end
