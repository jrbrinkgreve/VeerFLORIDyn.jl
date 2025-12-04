


# custom RP data generator for testing
function generate_RP_data_XY(nRPx::Int, nRPy::Int, par::Params, zloc::Float64)
    # define grid extents (adjust as you like)
    x_min, x_max = 0.0, 10.0 * par.D
    y_min, y_max = -1.5 * par.D, 1.5 * par.D

    xs = range(x_min, x_max; length = nRPx)
    ys = range(y_min, y_max; length = nRPy)

    nRP = nRPx * nRPy
    RP_coords = zeros(nRP, 3)

    k = 1
    for x in xs, y in ys
        RP_coords[k, 1] = x      # x coordinate
        RP_coords[k, 2] = y      # y coordinate
        RP_coords[k, 3] = zloc   # z coordinate (constant)
        k += 1
    end

    return nRP, RP_coords
end




# custom RP data generator for testing
function generate_RP_data_YZ(nRPx::Int, nRPy::Int, par::Params, xloc::Float64)
    # define grid extents (adjust as you like)
    y_min, y_max = -1.5 * par.D, 1.5 * par.D
    z_min, z_max = 0.0, 2.0 * par.D

    ys = range(y_min, y_max; length = nRPy)
    zs = range(z_min, z_max; length = nRPy)

    nRP = nRPx * nRPy
    RP_coords = zeros(nRP, 3)
    RP_coords[:,1] .= xloc   # x coordinate (constant)

    k = 1
    for y in ys, z in zs
        RP_coords[k, 2] = y      # y coordinate
        RP_coords[k, 3] = z   # z coordinate
        k += 1
    end

    return nRP, RP_coords
end


