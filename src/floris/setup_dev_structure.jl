#note: this file should be run after the exfiltrate operation



set.enable_veer = true
safehouse.set.enable_veer = true


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




function FLORISBuffers(n_pts::Int)
    return FLORISBuffers(
        Matrix{Float64}(undef, n_pts, 3),  # tmp_RPs
        Matrix{Float64}(undef, n_pts, 3),  # rotor_pts
        Vector{Float64}(undef, n_pts),     # sig_y
        Vector{Float64}(undef, n_pts),     # sig_z
        Vector{Float64}(undef, n_pts),     # x_0
        Matrix{Float64}(undef, n_pts, 2),  # delta
        Vector{Float64}(undef, n_pts),     # pc_y
        Vector{Float64}(undef, n_pts),     # pc_z
        Vector{Float64}(undef, n_pts),     # cw_y
        Vector{Float64}(undef, n_pts),     # cw_z
        Vector{Float64}(undef, n_pts),     # phi_cw
        Vector{Float64}(undef, n_pts),     # r_cw
        Vector{Bool}(undef, n_pts),        # core
        Vector{Bool}(undef, n_pts),        # nw
        Vector{Bool}(undef, n_pts),        # fw
        Vector{Float64}(undef, n_pts),     # tmp_RPs_r
        Vector{Float64}(undef, n_pts),     # gaussAbs
        Vector{Float64}(undef, n_pts),     # gaussWght
        Vector{Float64}(undef, n_pts),     # exp_y
        Vector{Float64}(undef, n_pts),     # exp_z
        Vector{Bool}(undef, n_pts),        # not_core

        n_pts,
        Matrix{Float64}(undef, n_pts, 3),  #rps_coords    
        #for getVeerVars! outputs
        Vector{Float64}(undef, n_pts),   #alpha, height dependent veer angle
        Vector{Float64}(undef, n_pts),   #beta
        Vector{Float64}(undef, n_pts),   #gamma
        Vector{Float64}(undef, n_pts),   #a_star
        Vector{Float64}(undef, n_pts),   #xi_0_hat
        Matrix{Float64}(undef, 3, 3),  #rotmtx
        Matrix{Float64}(undef, n_pts, 3),#coords_veered
        Vector{Float64}(undef, n_pts),   #shear_modifier
        Vector{Float64}(undef, n_pts),   #u_in_z
        Vector{Float64}(undef, n_pts),   #t_hat
        Vector{Float64}(undef, n_pts),   #sgn_t_hat
        Vector{Float64}(undef, n_pts),   #abs_t_hat
        Vector{Float64}(undef, n_pts),   #y_hat_c
        Vector{Float64}(undef, n_pts),   #y_c
        Vector{Float64}(undef, n_pts),   #theta
        Vector{Float64}(undef, n_pts),   #xi_0
        Vector{Float64}(undef, n_pts),   #xi_hat
        Vector{Float64}(undef, n_pts),   #chi
        Vector{Float64}(undef, n_pts),   #a
        Vector{Float64}(undef, n_pts),   #c1
        Vector{Float64}(undef, n_pts),   #c2
        Vector{Float64}(undef, n_pts),   #c3
        Vector{Float64}(undef, n_pts),   #c4
        Vector{Float64}(undef, n_pts),   #c5
        Vector{Float64}(undef, n_pts),   #c6
        Vector{Float64}(undef, n_pts),   #c7
        Vector{Float64}(undef, n_pts),   #xi
        Vector{Float64}(undef, n_pts),   #sigma
        Vector{Float64}(undef, n_pts),   #sigma_hat_squared
        Vector{Float64}(undef, n_pts),   #c
        Vector{Float64}(undef, n_pts),   #du
        Vector{Float64}(undef, n_pts),   #u
        Vector{Float64}(undef, n_pts),    #tmp

        Float64[],                         # T_red_arr (size set per call)
        Float64[],                         # T_aTI_arr (size set per call)
        Float64[],                         # T_Ueff (size 0 or 1)
        Float64[],                         # T_weight (size set per call)
    )
end



