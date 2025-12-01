using Plots
using ScatteredInterpolation


function plot_rotor_points(buf::bufferstruct, par::Params)
    x = buf.rps_coords[:, 1]
    y = buf.rps_coords[:, 2]

    scatter(x, y; scatter = true,
         xlabel = "x (m)",
         ylabel = "y (m)",
         title = "Rotor Points",
         legend = false)
end



function plot_velocity(buf::bufferstruct, par::Params)
    x = buf.rps_coords[:, 1]
    y = buf.rps_coords[:, 2]
    z = buf.u

    scatter(x, y;
            xlabel = "x (m)",
            ylabel = "y (m)",
            title = "Wake Velocity Field",
            markersize = 1,
            marker_z = z,
            color = :viridis,
            colorbar = true,
            legend = false)
end

function plot_velocity_advanced(buf::bufferstruct, par::Params)
    # 1. Prepare Data and Normalization
    D = 2 * par.R
    
    # Coordinates: x (streamwise, col 1) and y (lateral, col 2)
    x = buf.rps_coords[:, 1]
    y = buf.rps_coords[:, 2]

    # Data: Normalized Velocity Deficit (Delta U / U_h)
    # We assume buf.du is the velocity deficit (u_in - u)
    normalized_vel = buf.u ./ par.u_hub
    
    # 2. Build Interpolant
    pts = [x y]'                         # 2 × N (X and Y coordinates)
    itp = interpolate(Multiquadratic(), pts, normalized_vel)                
    # 3. Define Regular Grid for the X-Y plane
    nx, ny = 200, 200                     # adjust resolution
    xg = range(minimum(x), maximum(x), length = nx)
    yg = range(minimum(y), maximum(y), length = ny)

    X_grid = repeat(xg, ny)[:]
    Y_grid = repeat(yg', nx)[:]
    gridPts = [X_grid Y_grid]'

    # Evaluate the interpolant on the grid
    zg_interp = evaluate(itp, gridPts)
    
    # Reshape into a Matrix (X runs horizontally, Y runs vertically)
    ZG_MATRIX = reshape(zg_interp, ny, nx)
    ZG_MATRIX = ZG_MATRIX' # Transpose to ensure correct axis alignment

    # 4. Construct the Title String
    # Assume the horizontal slice is at hub height z_hub
    z_over_D = round(par.z_hub / D, digits=2)
    veer_val = par.alpha_gradient
    yaw_val  = par.beta

    title_str = "Normalized velocity U/u_h (Horizontal Slice at z/D ≈ $(z_over_D))\n" * "Veer: $(veer_val) °/m | Yaw: $(yaw_val)°"
    
    # 5. Plotting Configuration
    xg_norm = xg ./ D
    yg_norm = yg ./ D
    max_vel = maximum(normalized_vel)

    default(size = (800, 400), titlefontsize = 10)

    heatmap(
        xg_norm, yg_norm, ZG_MATRIX;
        xlabel = "x/D",           # Streamwise X is horizontal (corrected from your original code)
        ylabel = "y/D",           # Lateral Y is vertical (corrected from your original code)
        title  = title_str,
        colorbar_title = "U / U_h",
        c = :inferno,                    
        aspect_ratio = :equal,        
        clims = (0.0, max_vel * 1.05)
    )

end



function plot_velocity_advanced_YZ(buf::bufferstruct, par::Params  )


    #data
    x = buf.rps_coords[:, 2]
    y = buf.rps_coords[:, 3]
    z = buf.u

    # 1) build interpolant from scattered points
    pts = [x y]'                        # 2 × N
    itp = interpolate(Multiquadratic(), pts, z)

    # 2) define regular grid
    nx, ny = 200, 200             # adjust resolution
    xg = range(minimum(x), maximum(x), length = nx)
    yg = range(minimum(y), maximum(y), length = ny)

    X = repeat(xg, ny)[:]
    Y = repeat(yg', nx)[:]
    gridPts = [X Y]'

    zg = evaluate(itp, gridPts)
    ZG = reshape(zg, ny, nx)            # matrix for heatmap
    ZG = ZG'    #transpose to get axes right


    default(size = (700, 700))

    heatmap(
        xg, yg, ZG;
        xlabel = "y [m]",
        ylabel = "z [m]",
        title  = "Effective Wind Speed",
        colorbar_title = "WindSpeed [m s⁻¹]",
        c = :inferno,                    # or :viridis, etc."
        aspect_ratio = :10/3,

    )

end


using Plots, Printf

function plot_contour_YZ(buf::bufferstruct, par::Params)
    # 1. Define Diameter for Normalization
    D = par.R * 2.0

    # 2. Extract Data
    y_flat = buf.rps_coords[:, 2]
    z_flat = buf.rps_coords[:, 3]
    deficit_flat = buf.du ./ par.u_hub

    # 3. Extract Simulation Parameters for Title
    # We assume the slice is at a constant X, so we take the first point's X value
    x_val = buf.rps_coords[1, 1]
    x_over_D = round(x_val / D, digits=1)
    
    # Extract Veer and Yaw from Params (assuming they are in degrees based on your code)
    veer_val = par.alpha_gradient
    yaw_val  = par.beta

    # 4. Construct the Title String
    # \n creates a new line. We use @sprintf or string interpolation for formatting.
    # Note: If Veer/Yaw are floats, you might want to round them.
    title_str = "Normalized velocity deficit ΔU/u_h at x/D = $(x_over_D)\nVeer: $(veer_val) °/m | Yaw: $(yaw_val)°"

    # 5. Reconstruct Grid (Same as before)
    ys = sort(unique(y_flat))
    zs = sort(unique(z_flat))
    ny = length(ys)
    nz = length(zs)

    if length(deficit_flat) != ny * nz
        error("Buffer data size mismatch.")
    end

    grid_deficit = reshape(deficit_flat, ny, nz)

    # 6. Plot with New Title
    contourf(ys ./ D, zs ./ D, grid_deficit,
        xlabel = "y/D",
        ylabel = "z/D",
        title = title_str,        # <--- Updated Title
        titlefontsize = 10,       # Adjust font size if the two lines are too big
        fill = true,
        levels = 10,
        c = :inferno,
        aspect_ratio = :equal,
        clims = (0.0, maximum(deficit_flat))
    )
end


function plotting_surf(buf::bufferstruct, par::Params  )
    x = buf.rps_coords[:, 1]
    y = buf.rps_coords[:, 2]
    z = buf.u

    surface(x, y, z;
         xlabel = "x (m)",
         ylabel = "y (m)",
         zlabel = "u (m/s)",
         title = "Wake Velocity Field Surface",
         legend = false)


end




