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


#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------




function handle_single_turbine_veer!(buffers, RPl, RPw, location_t, set, windshear, d_rotor)
    # Avoid allocating RPl[:,3] and the broadcasted division by using a buffer

    nRP_local = size(RPl, 1)
    if length(buffers.tmp_RPs_r) < nRP_local
        error("FLORISBuffers.tmp_RPs_r too small: expected at least $(nRP_local) elements, got $(length(buffers.tmp_RPs_r)).\n" *
              "Ensure create_unified_buffers(.., floris) used the same rotor discretization.")
    end
    if set.shear_mode isa Shear_PowerLaw
        # Power law expects z normalized by hub height; clamp to > 0
        @inbounds for i in 1:nRP_local
            val = RPl[i, 3] / location_t[end, 3]
            buffers.tmp_RPs_r[i] = val > eps() ? val : eps()
        end
    else
        # Interpolation expects absolute height in meters
        @inbounds for i in 1:nRP_local
            buffers.tmp_RPs_r[i] = RPl[i, 3]
        end
    end
    z_view = @view buffers.tmp_RPs_r[1:nRP_local]
    redShear = getWindShearT(set.shear_mode, windshear, z_view)  #normalized!
    @infiltrate
    
    
    
    # Avoid allocating a view for RPw in the dot product
    acc = 0.0
    @inbounds for i in 1:nRP_local
        acc = muladd(RPw[i], redShear[i], acc)
    end
    T_red_scalar = acc
    # Persist result into buffers as 1-length arrays (optional use by callers)
    resize!(buffers.T_red_arr, 1); buffers.T_red_arr[1] = T_red_scalar
    resize!(buffers.T_aTI_arr, 0)
    resize!(buffers.T_Ueff, 0)
    resize!(buffers.T_weight, 0)
    return nothing

end



#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------




function setup_computation_buffers_veer!(buffers, nRP::Int, nT::Int)
    # Initialize outputs in buffers
    resize!(buffers.T_red_arr, nT); fill!(buffers.T_red_arr, 1.0)
    resize!(buffers.T_aTI_arr, max(nT - 1, 0))
    if nT > 1
        fill!(buffers.T_aTI_arr, 0.0)
    end
    resize!(buffers.T_weight, max(nT - 1, 0))
    if nT > 1
        fill!(buffers.T_weight, 0.0)
    end

    # Ensure buffers are properly sized
    if size(buffers.tmp_RPs, 1) < nRP
        error("Buffer tmp_RPs is too small: expected at least $(nRP) rows, got $(size(buffers.tmp_RPs, 1))")
    end
    if length(buffers.cw_y) < nRP
        error("Buffer arrays are too small: expected at least $(nRP) elements, got $(length(buffers.cw_y))")
    end
    
    # Use views of pre-allocated buffers to match the current discretization size exactly
    tmp_RPs = view(buffers.tmp_RPs, 1:nRP, :)
    sig_y = view(buffers.sig_y, 1:nRP)
    sig_z = view(buffers.sig_z, 1:nRP)
    x_0   = view(buffers.x_0, 1:nRP)
    delta = view(buffers.delta, 1:nRP, :)
    pc_y  = view(buffers.pc_y, 1:nRP)
    pc_z  = view(buffers.pc_z, 1:nRP)
    cw_y = view(buffers.cw_y, 1:nRP)
    cw_z = view(buffers.cw_z, 1:nRP)
    phi_cw = view(buffers.phi_cw, 1:nRP)
    r_cw = view(buffers.r_cw, 1:nRP)
    core = view(buffers.core, 1:nRP)
    nw = view(buffers.nw, 1:nRP)
    fw = view(buffers.fw, 1:nRP)
    tmp_RPs_r = view(buffers.tmp_RPs_r, 1:nRP)
    gaussAbs = view(buffers.gaussAbs, 1:nRP)
    gaussWght = view(buffers.gaussWght, 1:nRP)
    exp_y = view(buffers.exp_y, 1:nRP)
    exp_z = view(buffers.exp_z, 1:nRP)
    not_core = view(buffers.not_core, 1:nRP)
    
    
    
    #veer:
    rps_coords = view(buffers.rps_coords, :,:)
    alpha = view(buffers.alpha, :)
    beta = view(buffers.beta, :)
    gamma = view(buffers.gamma, :)
    a_star = view(buffers.a_star, :)
    xi_0_hat = view(buffers.xi_0_hat, :)
    rotmtx = view(buffers.rotmtx, :, :)
    coords_veered = view(buffers.coords_veered, :, :)
    shear_modifier = view(buffers.shear_modifier, :)
    u_in_z = view(buffers.u_in_z, :)
    t_hat = view(buffers.t_hat, :)
    sgn_t_hat = view(buffers.sgn_t_hat, :)
    abs_t_hat = view(buffers.abs_t_hat, :)
    y_hat_c = view(buffers.y_hat_c, :)
    y_c = view(buffers.y_c, :)
    theta = view(buffers.theta, :)
    xi_0 = view(buffers.xi_0, :)
    xi_hat = view(buffers.xi_hat, :)
    chi = view(buffers.chi, :)  
    a = view(buffers.a, :)
    c1 = view(buffers.c1, :)
    c2 = view(buffers.c2, :)
    c3 = view(buffers.c3, :)
    c4 = view(buffers.c4, :)
    c5 = view(buffers.c5, :)
    c6 = view(buffers.c6, :)
    c7 = view(buffers.c7, :)
    xi = view(buffers.xi, :)
    sigma = view(buffers.sigma, :)
    sigma_hat_squared = view(buffers.sigma_hat_squared, :)
    c = view(buffers.c, :)
    du = view(buffers.du, :)
    u = view(buffers.u, :)
    tmp = view(buffers.tmp, :)
    

    return (tmp_RPs, sig_y, sig_z, x_0, delta, pc_y, pc_z, cw_y, cw_z, phi_cw, r_cw, 
            core, nw, fw, tmp_RPs_r, gaussAbs, gaussWght, exp_y, exp_z, not_core, rps_coords,
            alpha, beta, gamma, a_star, xi_0_hat, rotmtx, coords_veered, shear_modifier, u_in_z,
            t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c, theta, xi_0, xi_hat, chi, a, c1, c2, c3, 
            c4, c5, c6, c7, xi, sigma, sigma_hat_squared, c, du, u, tmp)


end


#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------

function compute_wake_effects_veer!(buffers, views, iT, RPl, RPw, location_t, states_wf, 
                                states_t, d_rotor, floris, nRP, set, windshear)



    #get views:
    tmp_RPs, sig_y, sig_z, x_0, delta, pc_y, pc_z, cw_y, cw_z, phi_cw, r_cw, 
            core, nw, fw, tmp_RPs_r, gaussAbs, gaussWght, exp_y, exp_z, not_core, rps_coords,
            alpha, beta, gamma, a_star, xi_0_hat, rotmtx, coords_veered, shear_modifier, u_in_z,
            t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c, theta, xi_0, xi_hat, chi, a, c1, c2, c3, 
            c4, c5, c6, c7, xi, sigma, sigma_hat_squared, c, du, u, tmp = views
    


    #get yaw angle
    tmp_phi = size(states_wf,2) == 4 ? angSOWFA2world(states_wf[iT, 4]) :
                                       angSOWFA2world(states_wf[iT, 2])




    #construct relative RP locations: relevant for wake model evaluations
    for i in 1:nRP
        for j in 1:2   
            tmp_RPs[i, j] = RPl[i, j] - location_t[iT, j]   #relative location wrt to this turbine
        end         #tmp_RPs contains relative X,Y coordinates w.r.t nacelle of current turbine
    end             #as wake model uses absolute heights, Z coord is NOT relative    

    tmp_RPs[:, 3] = RPl[:, 3]   #absolute heights, no relative here

  

    
    
    
    #For this veered wake model, we use absolute z coordinates, wake model by
    #Mohammadi et al. uses this
    #as ground effects are considered there, we need to know absolute heights
    z_hub = location_t[iT, 3] #added for clarity
    #precompute for yaw frame alignment 
    cos_phi = cos(tmp_phi)    #precompute
    sin_phi = sin(tmp_phi)
    
    # Apply rotation matrix manually to avoid allocation
    #rotation matrix is for aligning with yaw 
    for i in 1:nRP
        x = tmp_RPs[i, 1]
        y = tmp_RPs[i, 2]
        z = tmp_RPs[i, 3]
        tmp_RPs[i, 1] = cos_phi * x + sin_phi * y
        tmp_RPs[i, 2] = -sin_phi * x + cos_phi * y
        tmp_RPs[i, 3] = z
    end

    #prob a break on a failed sanity check
    if tmp_RPs[1, 1] <= 10
        print("something went wrong, veer.jl line 220")
        return nothing
    end


    #load variability in readable names
    a_val = states_t[iT, 1]
    yaw_deg = states_t[iT, 2]
    yaw = -deg2rad(yaw_deg)
    TI = states_t[iT, 3]
    Ct = calcCt(a_val, yaw_deg)
    TI0 = states_wf[iT, 3]
    u_hub = states_wf[iT, 1]





    # Compute mean_x now, before tmp_RPs is reused for other data
    mean_x = 0.0
    @inbounds for i in 1:nRP
        mean_x += tmp_RPs[i, 1]
    end
    mean_x /= nRP

   
    
    #get Mohammadi model parameters / velocity field writting inplace to u
    getVars_veer!(  
        tmp_RPs,
        alpha, yaw, gamma, a_star, xi_0_hat, coords_veered, shear_modifier, u_in_z,
        t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c, theta, xi_0, xi_hat, chi, a, c1, c2, c3, 
        c4, c5, c6, c7, xi, sigma, sigma_hat_squared, c, du, u, tmp,
        Ct, TI, TI0, floris, d_rotor[iT],  #get rotor diameter of turbine causing the wake
        z_hub, set, windshear, u_hub[iT] )  #get hub velocity of turbine causing the wake
    

end





#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------



function getVars_veer!(rps_coords,
            alpha, beta, gamma, a_star, xi_0_hat, coords_veered, shear_modifier, u_in_z,
            t_hat, sgn_t_hat, abs_t_hat, y_hat_c, y_c, theta, xi_0, xi_hat, chi, a, c1, c2, c3, 
            c4, c5, c6, c7, xi, sigma, sigma_hat_squared, c, du, u, tmp,
            CT, ti, ti0, floris::Floris, D, z_hub, set, windshear, u_hub)

    # Parameters and constraints
    sqrt3 = sqrt(3.0)
    pisq = pi^2.0    
    pim1 = pi - 1.0
    nRP, _ = size(rps_coords)
    R = D / 2.0
  




    for i in 1:nRP                #use inbounds for performance
        #do the wake model evaluation per RP here:
        # determining veered wind directon at RP height
        alpha[i] = deg2rad(floris.veer_gradient * (rps_coords[i,3] - z_hub))  #linear veer profile
        gamma[i] = deg2rad(beta) + alpha[i]
        
        #eq 3    
        a_star[i] = (1.0  + sqrt(1.0 - CT * cos(gamma[i])^2)    ) / (2.0 * sqrt(1.0 - CT * cos(gamma[i])^2))
        
        xi_0_hat[i] = R * sqrt(a_star[i])
        
        
        #rotate RP coordinates to veered frame at each height
        coords_veered[i,1] = rps_coords[i,1] * cos(alpha[i]) + rps_coords[i,2] * sin(alpha[i])
        coords_veered[i,2] = rps_coords[i,1] * -sin(alpha[i]) + rps_coords[i,2] * cos(alpha[i])
        coords_veered[i,3] = rps_coords[i,3] 

        #compute shear modifier at RP height
        shear_modifier[i] = getWindShearT(set.shear_mode, windshear, coords_veered[i,3]/z_hub)
        u_in_z[i] = u_hub * shear_modifier[i]   #note: check whether this shear modifier is correct
        

        #AAAAA      CONTINUE DEVELOPMENT FROM HERE: ADD USTAR TO FLORIS PARAMS STRUCT

        t_hat[i] = (
            -1.44 * u_in_z[i]/ floris.u_star * R / xi_0_hat[i] * CT  
            * cos(gamma[i])^2 * sin(gamma[i])  *
            (1.0 - exp( -0.35   * floris.u_star / (u_in_z[i]) *  coords_veered[i,1] / R )   )
        )
        
        y_hat_c[i] = ((pim1 *  abs(t_hat[i])^3 + 2.0sqrt3 * pisq * t_hat[i]^2 + 48.0pim1^2 * abs(t_hat[i] ) ) / 
                (2.0*pi*pim1 * t_hat[i]^2 + 4.0sqrt3 * pisq * abs(t_hat[i]) + 96.0 * pim1^2) * sign(t_hat[i]) - 
                (2.0 / pi) * t_hat[i]  / (((coords_veered[i,3] + z_hub  ) / xi_0_hat[i])^2 - 1.0)  

        )

        y_c[i] = y_hat_c[i] * xi_0_hat[i]
        
        if abs(coords_veered[i,3] - z_hub) < 1e-10
             theta[i] = 0.0
        else
            theta[i] = atan( (coords_veered[i,3] - z_hub) / (coords_veered[i,2] - y_c[i]) )
        end


        #eq 9: initial wake shape
        xi_0[i] = R*sqrt(a_star[i]) * abs(cos(gamma[i])) / sqrt(1.0 - (sin(gamma[i]) * sin(theta[i]) )^2   )
        
        if gamma[i] < floris.angle_tol
            chi[i] = 0.0
            xi_hat[i] = 1.0
        else

            chi[i] = 1.0 / (floris.tsr * sin(gamma[i]))
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
            (floris.k * floris.u_star / (u_in_z[i])  * coords_veered[i,1] + 0.4 * xi_0_hat[i]) *
            (floris.k * floris.u_star / (u_in_z[i])  * coords_veered[i,1] + 0.4 * xi_0_hat[i] * cos(gamma[i]))
        ) # == (kx + 0.4 xi_0_hat     )(         ) 
        sigma[i] = (floris.k * floris.u_star / (u_in_z[i])  * coords_veered[i,1] + 0.4 * xi[i])

        #eq 15: velocity deficit
        c[i] = 1.0 - sqrt(max(0.001,      1.0 - R^2 * CT * cos(gamma[i])^3 / (2.0 * sigma_hat_squared[i]) ))
        du[i] = c[i] * exp(- ((coords_veered[i,2] - y_c[i])^2 + (coords_veered[i,3] - z_hub)^2 )  / (2.0*sigma[i]^2)   ) * (u_in_z[i])
        u[i] = u_in_z[i] - du[i]

        
        
    end
    @infiltrate

end




#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------




function compute_final_wind_shear_veer!(buffers, RPl, RPw, location_t, set::Settings, 
                                  windshear, tmp_RPs_r, states_wf)

end












#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------





