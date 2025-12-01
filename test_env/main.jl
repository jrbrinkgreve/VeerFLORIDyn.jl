include("structs.jl")
include("funcs.jl")
include("plotting_scripts.jl")
include("data_generators.jl")

#testbench for veer wake model

nRPx = 51
nRPy = 51
buf = bufferstruct(nRPx * nRPy) 
par = Params()


par.alpha_gradient = 0.05   #set alpha gradient for veering test
par.beta = 10.0              #set yaw angle for veering test

zloc = par.z_hub    #set z location for XY plane
xloc = 10.0 * par.D  #set x location for YZ plane

RP_data = generate_RP_data_YZ(nRPx, nRPy , par, xloc)
runFUNCTIONS!(buf, par, RP_data)
plot_contour_YZ(buf, par)