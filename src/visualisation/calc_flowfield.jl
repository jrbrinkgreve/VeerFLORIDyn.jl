# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Base.Threads

"""
    ThreadBuffers

Thread-local buffers for parallel flow field computation.

# Fields
- `thread_buffers::Vector{WindFarm}`: Thread-local WindFarm objects for each thread
- `thread_unified_buffers::Vector{UnifiedBuffers}`: Thread-local UnifiedBuffers for computation
"""
struct ThreadBuffers
    thread_buffers::Vector{WindFarm}
    thread_unified_buffers::Vector{UnifiedBuffers}
end

# UnifiedBuffers struct and create_unified_buffers function are defined in floridyn_cl/structs.jl

"""
    create_thread_buffers(wf::WindFarm, nth::Int, floris::Floris) -> ThreadBuffers

Create thread-local buffers for parallel flow field computation with FLORIS parameters.

This function pre-allocates all necessary data structures for each thread to avoid
race conditions and memory allocations during the parallel computation loop.

# Arguments
- `wf::WindFarm`: Original wind farm object to use as template
- `nth::Int`: Number of threads to create buffers for
- `floris::Floris`: FLORIS parameters for creating proper FLORIS buffers

# Returns
- `ThreadBuffers`: Struct containing all thread-local buffers

# Performance Notes
- Each thread gets its own copy of the WindFarm structure
- Pre-allocates all arrays to minimize allocations during computation
- Sets up dependency structure for virtual turbines at grid points
"""
function create_thread_buffers(wf::WindFarm, nth::Int, floris::Floris)
    # Pre-allocate unified buffers for each thread (contains GP)
    thread_unified_buffers = Vector{UnifiedBuffers}(undef, nth)
    thread_buffers = Vector{WindFarm}(undef, nth)

    for tid in 1:nth
        # Create unified buffers with proper FLORIS parameters for this thread
        ub = create_unified_buffers(wf, floris)
        thread_unified_buffers[tid] = ub
        # Use the prebuilt GP from unified buffers as the thread's WindFarm buffer
        thread_buffers[tid] = ub.gp
    end

    return ThreadBuffers(thread_buffers, thread_unified_buffers)
end

"""
    create_thread_buffers(wf::WindFarm, nth::Int) -> ThreadBuffers

Create thread-local buffers for parallel flow field computation.

This function pre-allocates all necessary data structures for each thread to avoid
race conditions and memory allocations during the parallel computation loop.

# Arguments
- `wf::WindFarm`: Original wind farm object to use as template
- `nth::Int`: Number of threads to create buffers for

# Returns
- `ThreadBuffers`: Struct containing all thread-local buffers

# Performance Notes
- Each thread gets its own copy of the WindFarm structure
- Pre-allocates all arrays to minimize allocations during computation
- Sets up dependency structure for virtual turbines at grid points
"""
function create_thread_buffers(wf::WindFarm, nth::Int)
    # Pre-allocate unified buffers for each thread (contains GP)
    thread_unified_buffers = Vector{UnifiedBuffers}(undef, nth)
    thread_buffers = Vector{WindFarm}(undef, nth)

    for tid in 1:nth
        # Create unified buffers with default rotor discretization
        ub = create_unified_buffers(wf, 50)
        thread_unified_buffers[tid] = ub
        # Use the prebuilt GP from unified buffers as the thread's WindFarm buffer
        thread_buffers[tid] = ub.gp
    end

    return ThreadBuffers(thread_buffers, thread_unified_buffers)
end

"""
    update_thread_buffers!(buffers::ThreadBuffers, wf::WindFarm) -> Nothing

Update wind field states in all thread-local wind farm buffers without allocating memory.

This function efficiently updates the wind field states (`States_WF`) and optional interpolation 
coefficients (`C_Vel`, `C_Dir`) in all thread-local WindFarm objects to match the current 
wind conditions from the master WindFarm object. This is useful when wind conditions change 
during simulation and the thread buffers need to be synchronized.

# Arguments
- `buffers::ThreadBuffers`: Thread-local buffers containing WindFarm copies for each thread
- `wf::WindFarm`: Master WindFarm object with updated wind field states

# Performance Notes
- Uses in-place assignment (`.=`) to avoid memory allocations
- Only updates wind-related fields, preserving other thread-specific modifications
- Updates all threads' buffers to maintain consistency
- Handles optional interpolation coefficient matrices when present

# Fields Updated
- `States_WF`: Wind field states matrix (velocity, direction, turbulence intensity)
- `C_Vel`: Velocity interpolation coefficients (if present in WindFarm type)  
- `C_Dir`: Direction interpolation coefficients (if present in WindFarm type)

# Example
```julia
# Create thread buffers
buffers = create_thread_buffers(wf, nthreads())

# ... wind conditions change ...

# Update all thread buffers with new wind states (non-allocating)
update_thread_buffers!(buffers, wf)

# Continue with flow field computation using updated buffers
Z = getMeasurements(buffers, X, Y, nM, zh, wf, set, floris, wind)
```

# See Also
- [`create_thread_buffers`](@ref): Create initial thread-local buffers
- [`getMeasurements`](@ref): Parallel flow field computation using thread buffers
"""
function update_thread_buffers!(buffers::ThreadBuffers, wf::WindFarm)
    # Update wind field states in all thread-local WindFarm objects
    for tid in 1:length(buffers.thread_buffers)
        thread_wf = buffers.thread_buffers[tid]
        
        # Update wind field states matrix (non-allocating in-place assignment)
        # Note: Only update the original turbine portion, not the extra grid point
        original_size = size(wf.States_WF)
        thread_wf.States_WF[1:original_size[1], 1:original_size[2]] .= wf.States_WF
        
        # Update optional interpolation coefficient matrices if they exist
        if hasfield(typeof(wf), :C_Vel) && isdefined(wf, :C_Vel)
            if hasfield(typeof(thread_wf), :C_Vel) && isdefined(thread_wf, :C_Vel)
                # Only update the original turbine portion
                original_nT = wf.nT
                thread_wf.C_Vel[1:original_nT, :] .= wf.C_Vel
            end
        end
        
        if hasfield(typeof(wf), :C_Dir) && isdefined(wf, :C_Dir)
            if hasfield(typeof(thread_wf), :C_Dir) && isdefined(thread_wf, :C_Dir)
                # Only update the original turbine portion
                original_nT = wf.nT
                thread_wf.C_Dir[1:original_nT, :] .= wf.C_Dir
            end
        end
    end
    
    return nothing
end



"""
    getMeasurements(buffers::ThreadBuffers, mx::Matrix, my::Matrix, nM::Int, zh::Real,
                    wf::WindFarm, set::Settings, floris::Floris, wind::Wind) -> Array{Float64,3}

Calculate flow field measurements at specified grid points by treating them as virtual turbines.

This function computes flow field properties (velocity reduction, added turbulence, effective wind speed)
at grid points by creating virtual turbines at each location and running the FLORIS wake model.
Each grid point is treated as a turbine that depends on all real turbines in the wind farm,
allowing wake effects to be captured in the flow field visualization.

# Arguments
- `buffers::ThreadBuffers`: Pre-allocated thread-local buffers created with [`create_thread_buffers`](@ref);
    for Julia 1.12 use `create_thread_buffers(wf, nthreads() + 1, floris)`; for single-thread use `create_thread_buffers(wf, 1, floris)`
- `mx::Matrix`: X-coordinates of grid points (m)
- `my::Matrix`: Y-coordinates of grid points (m)  
- `nM::Int`: Number of measurements to compute (typically 3)
- `zh::Real`: Hub height for measurements (m)
- `wf::WindFarm`: Wind farm object containing turbine data. See: [`WindFarm`](@ref)
  - `wf.nT`: Number of real turbines
  - `wf.StartI`: Starting indices for turbine data
  - `wf.posBase`, `wf.posNac`: Turbine positions
  - `wf.States_*`: Turbine state matrices
- `set::Settings`: Settings object containing simulation parameters. See: [`Settings`](@ref)
- `floris::Floris`: FLORIS model parameters for wake calculations. See: [`Floris`](@ref)
- `wind::Wind`: Wind field configuration. See: [`Wind`](@ref)

# Returns
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions `(size(mx,1), size(mx,2), nM)`
  - `mz[:,:,1]`: Velocity reduction
  - `mz[:,:,2]`: Added turbulence intensity
  - `mz[:,:,3]`: Effective wind speed

# Algorithm
For each grid point:
1. Creates a temporary wind farm with all original turbines plus one virtual turbine at the grid point
2. Sets the virtual turbine to depend on all real turbines (to capture wake effects)
3. Runs the FLORIS simulation to compute wake-affected flow properties
4. Extracts the result for the virtual turbine position

# Performance Notes
- Multi-threaded implementation using `@threads` for parallel processing of grid points when more than one buffer is provided
- With a single buffer (`length(buffers.thread_buffers) == 1`), runs in a single-threaded loop
- Each grid point requires a full wind farm simulation, so computation time scales with grid size
- Uses thread-local buffers created by [`create_thread_buffers`](@ref) to avoid race conditions
- On Julia 1.12 create `nthreads() + 1` buffers to accommodate thread indexing

# Example
```julia
# Create a 10x10 grid from 0 to 1000m
x_range = 0:100:1000
y_range = 0:100:1000
mx = repeat(collect(x_range)', length(y_range), 1)
my = repeat(collect(y_range), 1, length(x_range))

# Calculate 3 measurements at 90m hub height (single-thread)
buffers = create_thread_buffers(wind_farm, 1, floris_model)
mz = getMeasurements(buffers, mx, my, 3, 90.0, wind_farm, settings, floris_model, wind_config)

# Extract effective wind speed field
wind_speed_field = mz[:, :, 3]
```

# See Also
- [`calcFlowField`](@ref): Higher-level function that uses this to create complete flow field data
- [`setUpTmpWFAndRun!`](@ref): Underlying simulation function used for each grid point
"""
function getMeasurements(buffers, mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind)
    size_mx = size(mx)
    mz = zeros(size_mx[1], size_mx[2], nM)
    
    if length(buffers.thread_buffers) == 1
        # Single-threaded loop when only one buffer is provided
        for iGP in 1:length(mx)
            # Use the single available buffer
            GP = buffers.thread_buffers[1]
            unified_buffers = buffers.thread_unified_buffers[1]

            xGP = mx[iGP]
            yGP = my[iGP]

            GP.posBase[end, 1] = xGP
            GP.posBase[end, 2] = yGP
            GP.posBase[end, 3] = 0.0

            GP.posNac[end, 1] = 0.0
            GP.posNac[end, 2] = 0.0
            GP.posNac[end, 3] = zh

            # Reset the grid point state
            for j in 1:size(GP.States_T, 2)
                GP.States_T[end, j] = 0.0
            end

            interpolateOPs!(unified_buffers, GP.intOPs, GP)
            
            setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)
            tmpM = unified_buffers.M_buffer
        

            @views gridPointResult = tmpM[end, :]

            # Map linear index to (row, col) for column-major arrays
            q, r = divrem(iGP - 1, size_mx[1])
            row = r + 1
            col = q + 1
            @views mz[row, col, 1:3] .= gridPointResult
        end
    else
        # Parallel loop using @threads
        @threads :static for iGP in 1:length(mx)
            # Get thread-local buffers
            tid = threadid()
            GP = buffers.thread_buffers[tid]
            unified_buffers = buffers.thread_unified_buffers[tid]
            
            xGP = mx[iGP]
            yGP = my[iGP]
            
            # Thread-safe updates: create copies to avoid modifying shared arrays
            GP.posBase[end, 1] = xGP
            GP.posBase[end, 2] = yGP
            GP.posBase[end, 3] = 0.0
            
            GP.posNac[end, 1] = 0.0
            GP.posNac[end, 2] = 0.0
            GP.posNac[end, 3] = zh
            
            # Reset the grid point state (thread-safe element-wise assignment)
            for j in 1:size(GP.States_T, 2)
                GP.States_T[end, j] = 0.0
            end
            
            # Recalculate interpolated OPs for the updated geometry (non-allocating)
            interpolateOPs!(unified_buffers, GP.intOPs, GP)

            setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)
            tmpM = unified_buffers.M_buffer
            
            # Extract only the result for the grid point (last "turbine")
            @views gridPointResult = tmpM[end, :]
            
            # Convert linear index to subscripts (thread-safe)
            # Map linear index to (row, col) for column-major arrays
            q, r = divrem(iGP - 1, size_mx[1])
            row = r + 1
            col = q + 1
            
            # Thread-safe assignment using @views to avoid race conditions
            @views mz[row, col, 1:3] .= gridPointResult
        end
    end

    return mz
end

"""
    calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris; plt=nothing)

Generate full flow field plot data by calculating measurements across a grid.

This function creates a rectangular grid over the wind farm domain and calculates flow field
properties at each grid point by treating them as virtual turbines. The computation can be
performed in parallel if `set.threading` is true.

# Arguments
- `set::Settings`: Settings object containing simulation parameters
  - `set.threading`: If true, uses multi-threaded computation with `@threads`
  - `set.parallel`: If true, enables parallel-specific optimizations
- `wf::WindFarm`: Wind farm object containing turbine data
- `wind::Wind`: Wind field configuration  
- `floris::Floris`: FLORIS model parameters

# Keyword Arguments
- `plt=nothing`: Plot object for garbage collection control. If provided and `set.parallel` is true,
  automatically calls `plt.GC.enable(false)` before multithreading and `plt.GC.enable(true)` 
  after completion to prevent PyCall-related segmentation faults during parallel execution with
  ControlPlots loaded. To take full advantage of multithreading, executed the plotting in a separate process.
- `vis=nothing`: Visualization configuration object containing field limits and resolution settings.
  If provided, uses `vis.field_limits_min`, `vis.field_limits_max`, and `vis.field_resolution` 
  to define the computational grid. If not provided, defaults to domain [0,0,0] to [3000,3000,400] 
  meters with 20m resolution.

# Returns
- `Z::Array{Float64,3}`: 3D array of flow field measurements with dimensions `(ny, nx, 3)`
  - `Z[:,:,1]`: Velocity reduction factor
  - `Z[:,:,2]`: Added turbulence intensity  
  - `Z[:,:,3]`: Effective wind speed (m/s)
- `X::Matrix{Float64}`: X-coordinate grid (m)
- `Y::Matrix{Float64}`: Y-coordinate grid (m)

# Notes
- Grid resolution and domain are configurable via the `vis` parameter, or use default values for backward compatibility
- Hub height is taken from the first turbine in the wind farm

# Example
```julia
# Calculate flow field with threading and GC control
set.threading = true
Z, X, Y = calcFlowField(set, wf, wind, floris; plt)

# Extract velocity reduction field
velocity_reduction = Z[:, :, 1]

# Extract effective wind speed field  
wind_speed = Z[:, :, 3]
```

# See Also
- [`getMeasurements`](@ref): Function used internally to compute the flow field
- [`plotFlowField`](@ref): Visualization function for the generated data
"""


function calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris;
                       plt=nothing, vis=nothing)
    # Preallocate field
    nM = 3

    # Use vis struct fields if provided, otherwise fall back to defaults
    if vis !== nothing
        fieldLims = [vis.field_limits_min';
                     vis.field_limits_max']  # [xmin ymin zmin; xmax ymax zmax]
        fieldRes = vis.field_resolution
    else
        # Default values for backward compatibility
        fieldLims = [0.0 0.0 0.0;
                     3000.0 3000.0 400.0]  # [xmin ymin zmin; xmax ymax zmax]
        fieldRes = 20.0  # Resolution of the field in m
    end
    
    xAx = fieldLims[1,1]:fieldRes:fieldLims[2,1]
    yAx = fieldLims[1,2]:fieldRes:fieldLims[2,2]

    # Create coordinate grids (Julia equivalent of meshgrid)
    X = repeat(collect(xAx)', length(yAx), 1)
    Y = repeat(collect(yAx), 1, length(xAx))

    # Get hub height from first turbine
    zh = wf.posNac[1, 3]  #AAAAAAA modify this to plot at different heights to look at veer effects
    #@infiltrate
    
    # Get data
    if set.threading
        # Disable garbage collection before multithreading if not parallel
        if ! set.parallel
            GC.enable(false)
        end
        try
            # Create thread-local buffers using the new function with FLORIS parameters
            buffers = create_thread_buffers(wf, nthreads() + 1, floris)
            Z = getMeasurements(buffers, X, Y, nM, zh, wf, set, floris, wind)
        finally
            # Re-enable garbage collection after multithreading if not parallel
            if ! set.parallel
                GC.enable(true)
            end
        end
    else
        # Single-threaded path using getMeasurements with a single buffer
        buffers = create_thread_buffers(wf, 1, floris)
        Z = getMeasurements(buffers, X, Y, nM, zh, wf, set, floris, wind)
    end

    return Z, X, Y
end

"""
    calcFlowFieldCrossSection(set::Settings, wf::WindFarm, wind::Wind, floris::Floris;
                              fixed::Real, vis=nothing, orientation::Symbol=:NS)

Generate flow field data for a cross-section at a fixed location.

# Arguments
- `set::Settings`: Simulation settings
- `wf::WindFarm`: Wind farm object
- `wind::Wind`: Wind field configuration
- `floris::Floris`: FLORIS model parameters

# Keyword Arguments
- `fixed::Real`: The fixed coordinate (m) at which to take the cross-section
- `vis=nothing`: Visualization config with `field_limits_min`, `field_limits_max`,
  `field_resolution`. If not provided, defaults to the free axis ∈[0,3000], z∈[0,400] at 20m resolution.
- `orientation::Symbol=:NS`: Orientation of the cross-section:
  - `:NS` → fix x, vary y (South-North) and z; bottom axis is South-North
  - `:WE` → fix y, vary x (West-East) and z; bottom axis is West-East

# Returns
- `Z::Array{Float64,3}`: Flow field with dimensions `(nz, na, 3)`
  - `Z[:,:,1]`: Velocity reduction
  - `Z[:,:,2]`: Added turbulence intensity
  - `Z[:,:,3]`: Effective wind speed (m/s)
- `A::Matrix{Float64}`: Horizontal coordinate grid (y if `:NS`, x if `:WE`) (m)
- `Zh::Matrix{Float64}`: Z-coordinate grid (m)

# Example
```julia
# South-North cross-section at x=1000m
Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=1000.0, orientation=:NS)

# West-East cross-section at y=1500m
Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=1500.0, orientation=:WE)
```

# See Also
- [`plotFlowFieldCrossSection`](@ref): Plotting function for the generated data
- [`calcFlowField`](@ref): XY flow field equivalent
"""
function calcFlowFieldCrossSection(set::Settings, wf::WindFarm, wind::Wind, floris::Floris;
                                   fixed::Real, vis=nothing, xlims=nothing, orientation::Symbol=:NS)
    @assert orientation in (:NS, :WE) "orientation must be :NS or :WE"
    nM = 3

    if vis !== nothing
        a_min = orientation == :NS ? vis.field_limits_min[2] : vis.field_limits_min[1]
        a_max = orientation == :NS ? vis.field_limits_max[2] : vis.field_limits_max[1]
        z_min = vis.field_limits_min[3]
        z_max = vis.field_limits_max[3]
        fieldRes = vis.field_resolution
    else
        a_min, a_max = xlims !== nothing ? (xlims[1], xlims[2]) : (0.0, 3000.0)
        z_min, z_max = 0.0, 400.0
        fieldRes = 20.0
    end

    aAx = a_min:fieldRes:a_max
    zAx = z_min:fieldRes:z_max

    # Grid: rows = z-levels, cols = horizontal positions
    A  = repeat(collect(aAx)', length(zAx), 1)
    Zh = repeat(collect(zAx),  1, length(aAx))

    # Single buffer — no threading
    buffers        = create_thread_buffers(wf, 1, floris)
    GP             = buffers.thread_buffers[1]
    unified_buffers = buffers.thread_unified_buffers[1]

    size_grid = size(A)
    Z = zeros(size_grid[1], size_grid[2], nM)

    for iGP in 1:length(A)
        aGP = A[iGP]
        zGP = Zh[iGP]

        if orientation == :NS       # fix x, vary y
            GP.posBase[end, 1] = fixed
            GP.posBase[end, 2] = aGP
        else                        # fix y, vary x
            GP.posBase[end, 1] = aGP
            GP.posBase[end, 2] = fixed
        end
        GP.posBase[end, 3] = 0.0
        GP.posNac[end, 1]  = 0.0
        GP.posNac[end, 2]  = 0.0
        GP.posNac[end, 3]  = zGP

        for j in 1:size(GP.States_T, 2)
            GP.States_T[end, j] = 0.0
        end

        interpolateOPs!(unified_buffers, GP.intOPs, GP)
        setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)

        tmpM = unified_buffers.M_buffer
        @views gridPointResult = tmpM[end, :]

        q, r  = divrem(iGP - 1, size_grid[1])
        row   = r + 1
        col   = q + 1
        @views Z[row, col, 1:3] .= gridPointResult
    end

    return Z, A, Zh
end

"""
    plotFlowFieldCrossSection(plt, wf, A, Zh, Z, vis, fixed::Real, t=nothing;
                              msr::MSR=VelReduction, fig=nothing, orientation::Symbol=:NS)

Plot a cross-section of the flow field, mirroring the style of [`plotFlowField`](@ref).

# Arguments
- `plt`: Plotting package (e.g., ControlPlots.plt)
- `wf`: Wind farm object
- `A::Matrix`: Horizontal coordinate grid (y if `:NS`, x if `:WE`) (m)
- `Zh::Matrix`: Z-coordinate grid (m)
- `Z::Array{Float64,3}`: Flow field from [`calcFlowFieldCrossSection`](@ref), dims `(nz, na, 3)`
- `vis::Vis`: Visualization settings
- `fixed::Real`: The fixed coordinate of the cross-section, used in the title
- `t`: Optional time value for the plot title
- `msr::MSR`: Which measurement to plot. See: [`MSR`](@ref)
- `fig`: Optional figure name override
- `orientation::Symbol=:NS`: Must match the orientation used in [`calcFlowFieldCrossSection`](@ref):
  - `:NS` → bottom axis is South-North
  - `:WE` → bottom axis is West-East

# Returns
- `nothing`

# Example
```julia
# South-North cross-section at x=1000m
Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=1000.0, orientation=:NS)
plotFlowFieldCrossSection(plt, wf, A, Zh, Z, vis, 1000.0; orientation=:NS)

# West-East cross-section at y=1500m
Z, A, Zh = calcFlowFieldCrossSection(set, wf, wind, floris; fixed=1500.0, orientation=:WE)
plotFlowFieldCrossSection(plt, wf, A, Zh, Z, vis, 1500.0; orientation=:WE)
```

# See Also
- [`calcFlowFieldCrossSection`](@ref): Computes the cross-section data
- [`plotFlowField`](@ref): XY equivalent
"""
function plotFlowFieldCrossSection(plt, wf, A, Zh, Z, vis, fixed::Real, t=nothing;
                                   msr::MSR=VelReduction, fig=nothing, orientation::Symbol=:NS)
    @assert orientation in (:NS, :WE) "orientation must be :NS or :WE"

    mz_2d = Z[:, :, Int(msr)]

    if msr == VelReduction
        figure_name = "Velocity Reduction (Cross-Section)"
        label       = "Relative Wind Speed [%]"
        mz_2d .*= 100
        lev_min = vis.rel_v_min; lev_max = vis.rel_v_max
    elseif msr == AddedTurbulence
        figure_name = "Added Turbulence (Cross-Section)"
        label       = "Added Turbulence [%]"
        mz_2d .*= 100
        lev_min = 0.0; lev_max = vis.turb_max
    elseif msr == EffWind
        figure_name = "Effective Wind Speed (Cross-Section)"
        label       = "Wind speed [m/s]"
        lev_min = vis.v_min; lev_max = vis.v_max
    end
    if !isnothing(fig)
        figure_name = fig
    end

    # Axis labels and title depend on orientation
    if orientation == :NS
        xlabel_str = "South-North [m]"
        fixed_axis = "x"
    else
        xlabel_str = "West-East [m]"
        fixed_axis = "y"
    end

    title = figure_name * ", $(fixed_axis) = $(round(Int, fixed)) m"
    if t !== nothing
        time_str = lpad(string(round(Int, t)), 4, '0')
        title = title * ", t: " * time_str * " s"
    end

    aVec = A[1, :]    # unique horizontal values (one per column)
    zVec = Zh[:, 1]   # unique z values (one per row)

    fig_obj = plt.figure(figure_name, figsize=(7.25 * 0.84, 6 * 0.84))
    ax = plt.gca()
    n = 40
    levels = range(lev_min, stop=lev_max, length=n+1)

    plt.contourf(aVec, zVec, mz_2d, n; levels=levels, cmap="inferno")
    cb = plt.colorbar()
    cb.set_label(label, labelpad=3)

    plt.xlim(minimum(aVec), maximum(aVec))
    plt.ylim(minimum(zVec), maximum(zVec))
    plt.xlabel(xlabel_str)
    plt.ylabel("Height [m]")
    plt.title(title)
    plt.tight_layout(pad=1.0, h_pad=0.3, w_pad=0.3, rect=[0, 0, 1, 0.97])

    if vis.show_plots
        plt.show(block=false)
    end

    # Save if requested
    if vis.save && !vis.unit_test
        directory = vis.online ? vis.video_path : vis.output_path
        msr_name  = msr == VelReduction   ? "velocity_reduction" :
                    msr == AddedTurbulence ? "added_turbulence"   : "wind_speed"
        fixed_str = lpad(string(round(Int, fixed)), 4, '0')
        orient_str = orientation == :NS ? "x" : "y"
        filename  = t !== nothing ?
            joinpath(directory, "ff_cs_$(msr_name)_$(orient_str)$(fixed_str)m_t$(lpad(string(round(Int,t)),4,'0'))s.png") :
            joinpath(directory, "ff_cs_$(msr_name)_$(orient_str)$(fixed_str)m.png")
        try
            plt.savefig(filename, dpi=150, bbox_inches="tight", pad_inches=0.1, facecolor="white")
            if vis.print_filenames
                @info "Saving $filename"
            end
        catch e
            @debug "Failed to save cross-section plot" exception=(e, catch_backtrace())
        end
    end

    if vis.unit_test && vis.show_plots
        plt.pause(1.0)
        try plt.close(fig_obj) catch end
    end

    return nothing
end