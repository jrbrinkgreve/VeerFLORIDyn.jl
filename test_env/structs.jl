#structs for memory buffers

mutable struct teststruct
    alpha::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    a_star::Vector{Float64}
    xi_0_hat::Vector{Float64}
    rotmtx::Matrix{Float64}
    coords_veered::Matrix{Float64}
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

end


