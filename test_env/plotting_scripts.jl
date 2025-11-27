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


function plot_velocity_advanced(buf::bufferstruct, par::Params  )


    #data
    x = buf.rps_coords[:, 1]
    y = buf.rps_coords[:, 2]
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


    #default(fontfamily = "sans", size = (700, 700))

    heatmap(
        xg, yg, ZG;
        xlabel = "y [m]",
        ylabel = "x [m]",
        title  = "Effective Wind Speed",
        colorbar_title = "WindSpeed [m s⁻¹]",
        c = :jet,                    # or :viridis, etc."
        aspect_ratio = :10/3,

    )

end









