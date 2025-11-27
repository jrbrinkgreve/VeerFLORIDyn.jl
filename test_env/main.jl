include("structs.jl")
include("funcs.jl")
include("plotting_scripts.jl")

using Debugger

nRP = 50_00
buf = bufferstruct(nRP) 
par = Params()
RP_data = generate_RP_data(nRP, par)
runFUNCTIONS!(buf, par, RP_data)
plot_velocity_advanced(buf, par)

