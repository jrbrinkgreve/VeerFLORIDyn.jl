include("structs.jl")
include("funcs.jl")

using Debugger

nRP = 100
buf = bufferstruct(nRP) 
par = Params()
RP_data = generate_RP_data(nRP, par)


 
runFUNCTIONS!(buf, par, RP_data)



