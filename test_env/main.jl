include("structs.jl")
include("funcs.jl")


nRP = 100
buf = bufferstruct(nRP) 
par = Params()
RP_data = generate_RP_data(nRP)



runFUNCTIONS!(buf, par, RP_data)




