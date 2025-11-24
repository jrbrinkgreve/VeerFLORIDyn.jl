#funcs.jl



#constructor for bufferstruct
function bufferstruct(nRP::Int)::bufferstruct
    return bufferstruct(
        #input vars
        nRP,
        Matrix{Float64}(undef, 3, nRP),  #rps_coords    
        #for getVeerVars! outputs
        Vector{Float64}(undef, nRP),   #alpha, height dependent veer angle
        Vector{Float64}(undef, nRP),   #beta
        Vector{Float64}(undef, nRP),   #gamma
        Vector{Float64}(undef, nRP),   #a_star
        Vector{Float64}(undef, nRP),   #xi_0_hat
        Matrix{Float64}(undef, 3, 3),  #rotmtx
        Matrix{Float64}(undef, 3, nRP),#coords_veered
        Vector{Float64}(undef, nRP),   #t_hat
        Vector{Float64}(undef, nRP),   #sgn_t_hat
        Vector{Float64}(undef, nRP),   #abs_t_hat
        Vector{Float64}(undef, nRP),   #y_hat_c
        Vector{Float64}(undef, nRP),   #y_c
        Vector{Float64}(undef, nRP),   #theta
        Vector{Float64}(undef, nRP),   #xi_0
        Vector{Float64}(undef, nRP),   #xi_hat
        Vector{Float64}(undef, nRP),   #chi
        Vector{Float64}(undef, nRP),   #a
        Vector{Float64}(undef, nRP),   #c1
        Vector{Float64}(undef, nRP),   #c2
        Vector{Float64}(undef, nRP),   #c3
        Vector{Float64}(undef, nRP),   #c4
        Vector{Float64}(undef, nRP),   #c5
        Vector{Float64}(undef, nRP),   #c6
        Vector{Float64}(undef, nRP),   #c7
        Vector{Float64}(undef, nRP),   #xi
        Vector{Float64}(undef, nRP),   #sigma
        Vector{Float64}(undef, nRP),   #sigma_hat_squared
        Vector{Float64}(undef, nRP),   #c
        Vector{Float64}(undef, nRP),   #du
        Vector{Float64}(undef, nRP),   #u
        #output arrays here
            

    )
end





#constructor for Params struct
function Params()::Params
    return Params(
        0.05,
    )
end


#custom RP data generator for testing
function generate_RP_data(nRP::Int)
    RP_coords = 100*rand(3, nRP)  #example: random coordinates for rps_coords
    return nRP, RP_coords
end




#gets views: essentially pointers to the preallocated arrays in the bufferstruct    
function setup_computation_buffers!(buf::bufferstruct, nRP::Int)
    rps_coords = view(buf.rps_coords, :,:)
    alpha = view(buf.alpha, :)
    beta = view(buf.beta, :)
    gamma = view(buf.gamma, :)
    a_star = view(buf.a_star, :)
    xi_0_hat = view(buf.xi_0_hat, :)
    rotmtx = view(buf.rotmtx, :, :)
    coords_veered = view(buf.coords_veered, :, :)
    t_hat = view(buf.t_hat, :)
    sgn_t_hat = view(buf.sgn_t_hat, :)
    abs_t_hat = view(buf.abs_t_hat, :)
    y_hat_c = view(buf.y_hat_c, :)
    y_c = view(buf.y_c, :)
    theta = view(buf.theta, :)
    xi_0 = view(buf.xi_0, :)
    xi_hat = view(buf.xi_hat, :)
    chi = view(buf.chi, :)
    a = view(buf.a, :)
    c1 = view(buf.c1, :)
    c2 = view(buf.c2, :)
    c3 = view(buf.c3, :)
    c4 = view(buf.c4, :)
    c5 = view(buf.c5, :)
    c6 = view(buf.c6, :)
    c7 = view(buf.c7, :)
    xi = view(buf.xi, :)
    sigma = view(buf.sigma, :)
    sigma_hat_squared = view(buf.sigma_hat_squared, :)
    c = view(buf.c, :)
    du = view(buf.du, :)
    u = view(buf.u, :)
    return (nRP, rps_coords, alpha, beta, gamma, a_star, xi_0_hat, rotmtx,
            coords_veered, t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c,
            theta, xi_0, xi_hat, chi, a, c1, c2, c3, c4, c5, c6, c7,
            xi, sigma, sigma_hat_squared, c, du, u)
end




function compute_wake_effects!(buf::bufferstruct, par::Params, views, RP_data)
    nRP, rps_coords, alpha, beta, gamma, a_star, xi_0_hat, rotmtx,
    coords_veered, t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c,
    theta, xi_0, xi_hat, chi, a, c1, c2, c3, c4, c5, c6, c7,
    xi, sigma, sigma_hat_squared, c, du, u = views
   
    for i in 1:nRP
        #do the wake model evaluation per RP here
        rps_coords = RP_data[2]
        print(rps_coords)


    end

end











function runFUNCTIONS!(buf::bufferstruct, par::Params, RP_data)
    nRP = RP_data[1]
    views = setup_computation_buffers!(buf, nRP)
    compute_wake_effects!(buf, par, views, RP_data)
end



