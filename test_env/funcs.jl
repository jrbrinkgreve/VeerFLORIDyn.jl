#funcs.jl
using Debugger
using ControlPlots

#constructor for bufferstruct
function bufferstruct(nRP::Int)::bufferstruct
    return bufferstruct(
        #input vars
        nRP,
        Matrix{Float64}(undef, nRP, 3),  #rps_coords    
        #for getVeerVars! outputs
        Vector{Float64}(undef, nRP),   #alpha, height dependent veer angle
        Vector{Float64}(undef, nRP),   #beta
        Vector{Float64}(undef, nRP),   #gamma
        Vector{Float64}(undef, nRP),   #a_star
        Vector{Float64}(undef, nRP),   #xi_0_hat
        Matrix{Float64}(undef, 3, 3),  #rotmtx
        Matrix{Float64}(undef, nRP, 3),#coords_veered
        Vector{Float64}(undef, nRP),   #shear_modifier
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
        
            

    )
end




function Params()::Params
    return Params(      
        #example parameter values from Mohammadi et al.
        0.05,    #alpha_gradient
        63.0,   #R                          #import these later from turbine data
        126.0,  #D
        90.0,  #z_hub
        8.0,    #u_hub
        0.4,    #u_star 
        0.6,    #k
        0.66,    #CT
        1e-10,   #angle_tolerance for chi computation
        8.0,     #lambda, TSR
        15.0,     #beta, yaw angle in degrees
    )
end


#custom RP data generator for testing
function generate_RP_data(nRP::Int, par::Params)
    RP_coords = zeros(nRP, 3);   
    RP_coords[:,1] = 10.0* par.D * rand(nRP,1)           #example: random coordinates for rps_coords, tall matrix with 3 cols [x1, y1, z1]
    RP_coords[:,2] = 6.0 * par.D * (rand(nRP,1) .- 0.5)  # y deviation                                                        [x2, y2, z2] etc
    RP_coords[:,3] = 1.0 * par.D * (rand(nRP,1) .- 0.5)  # z deviation  
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
    shear_modifier = view(buf.shear_modifier, :)
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
            coords_veered, shear_modifier, t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c,
            theta, xi_0, xi_hat, chi, a, c1, c2, c3, c4, c5, c6, c7,
            xi, sigma, sigma_hat_squared, c, du, u)
end





@inline function windshearModifierPlaceholder(z::Float64)::Float64
    #placeholder for wind shear profile
    return  (1.0 + 0.001 * z) #0.1% per meter increase
end



function compute_wake_effects!(buf::bufferstruct, par::Params, views, RP_data)
    nRP, rps_coords, alpha, beta, gamma, a_star, xi_0_hat, rotmtx,
    coords_veered, shear_modifier, t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c,
    theta, xi_0, xi_hat, chi, a, c1, c2, c3, c4, c5, c6, c7,
    xi, sigma, sigma_hat_squared, c, du, u = views #use pointers to buffer
    rps_coords = RP_data[2]
    sqrt3 = sqrt(3.0)
    pisq = pi^2
    pim1 = pi - 1.0


    
    @inbounds for i in 1:nRP                                   #add @inbounds later!
        #do the wake model evaluation per RP here:
        
  

        # determining veered wind directon at RP height
        alpha[i] = deg2rad(par.alpha_gradient * rps_coords[i,3])  #linear veer profile
        gamma[i] = deg2rad(par.beta) + alpha[i]

        #eq 3    
        a_star[i] = (1.0  + sqrt(1.0 - par.CT * cos(gamma[i])^2)    ) / (2.0 * sqrt(1.0 - par.CT * cos(gamma[i])^2))
        xi_0_hat[i] = par.R * sqrt(a_star[i])
        
        
        #rotate RP coordinates to veered frame
        coords_veered[i,1] = rps_coords[i,1] * cos(alpha[i]) + rps_coords[i,2] * sin(alpha[i])
        coords_veered[i,2] = rps_coords[i,1] * -sin(alpha[i]) + rps_coords[i,2] * cos(alpha[i])
        coords_veered[i,3] = rps_coords[i,3] 
        
        #compute shear modifier at RP height
        shear_modifier[i] = windshearModifierPlaceholder(coords_veered[i,3])
        
        
        #compute t_hat
        t_hat[i] = (
            -1.44 * shear_modifier[i] * par.u_hub / par.u_star * par.R / xi_0_hat[i] * par.CT  
            * cos(gamma[i])^2 * sin(gamma[i])  *
            (1.0 - exp( -0.35   * par.u_star / (shear_modifier[i] * par.u_hub) *  coords_veered[i,1] / par.R )   )
        )
        y_hat_c[i] = ((pim1 *  abs(t_hat[i])^3 + 2.0sqrt3 * pisq * t_hat[i]^2 + 48.0pim1^2 * abs(t_hat[i] ) ) / 
                (2.0*pi*pim1 * t_hat[i]^2 + 4.0sqrt3 * pisq * abs(t_hat[i]) + 96.0 * pim1^2) * sign(t_hat[i]) - 
                (2.0 / pi) * t_hat[i]  / (((coords_veered[i,3] + 2.0*par.z_hub  ) / xi_0_hat[i])^2 - 1.0)  

        )     #note: ... + 2.0*par.z_hub  ) ... as z is defined from nacelle height, not ground level
       
        y_c[i] = y_hat_c[i] * xi_0_hat[i]

        theta[i] = atan( coords_veered[i,3] / (coords_veered[i,2] - y_c[i])  )



        #eq 9: initial wake shape
        xi_0[i] = par.R*sqrt(a_star[i]) * abs(cos(theta[i])) / sqrt(1.0 - (sin(gamma[i]) * sin(theta[i]) )^2   )
        
        if gamma[i] < par.angle_tolerance
            chi[i] = 0.0
            xi_hat[i] = 1.0
        else
            # Assume t_hat, gamma, par are defined and indexed

            chi[i] = 1.0 / (par.lambda * sin(gamma[i]))
            a[i] = 1.263 * cos(0.33 * chi[i])
            
            c1[i] = 0.5 * tanh(t_hat[i]^2 / (4.0*a[i]))
            c2[i] = (-1.0/3.0) * tanh(t_hat[i]^3 / (8.0*a[i]))
            c3[i] = -0.25 * tanh(t_hat[i]^3 / (8.0*a[i]))
            c4[i] = (-1.0/6.0) * tanh(t_hat[i]^4 / (16.0*a[i]))
            c5[i] = (5.0/16.0) * tanh(t_hat[i]^4 / (16.0*a[i]))
            c6[i] = (-5.0/48.0) * tanh(t_hat[i]^4 / (16.0*a[i]))
            c7[i] = (7.0/48.0) * tanh(t_hat[i]^4 / (16.0*a[i]))

            xi_hat[i] = 1.0 - a[i]*(
                c1[i] * cos(2.0theta[i]) +
                c2[i] * chi[i] * sin(2.0theta[i]) +
                c3[i] * cos(3.0theta[i]) +
                c4[i] * chi[i]^2 * cos(2.0theta[i]) +
                c5[i] * chi[i] * sin(3.0theta[i]) +
                c6[i] * cos(2.0theta[i]) +
                c7[i] * cos(4.0theta[i])
            )
        end
        xi[i] = xi_0[i] * xi_hat[i]
        sigma_hat_squared[i] = (
            (par.k * par.u_star / (par.u_hub * shear_modifier[i])  * coords_veered[i,1] + 0.4 * xi_0_hat[i]) *
            (par.k * par.u_star / (par.u_hub * shear_modifier[i])  * coords_veered[i,1] + 0.4 * xi_0_hat[i] * cos(gamma[i]))
        ) # == (kx + 0.4 xi_0_hat     )(         ) 
        sigma[i] = (par.k * par.u_star / (par.u_hub * shear_modifier[i])  * coords_veered[i,1] + 0.4 * xi[i])


        #eq 15: velocity deficit
        c[i] = 1.0 - sqrt(max(0.001,      1.0 - par.R^2 * par.CT * cos(gamma[i])^3 / (2.0 * sigma_hat_squared[i]) ))
        du[i] = c[i] * exp(- ((coords_veered[i,2] - y_c[i])^2 + coords_veered[i,3] )^2  / (2.0sigma[i]^2)   )
        u[i] = par.u_hub * shear_modifier[i] - du[i]
    end
    #print(size(du))

end








 


function runFUNCTIONS!(buf::bufferstruct, par::Params, RP_data)
    nRP = RP_data[1]
    views = setup_computation_buffers!(buf, nRP)
    compute_wake_effects!(buf, par, views, RP_data)
    plot
end



