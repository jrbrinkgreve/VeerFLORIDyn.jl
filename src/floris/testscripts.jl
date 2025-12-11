#=
function runFLORIS!(buffers::FLORISBuffers, set::Settings, location_t, states_wf, states_t, d_rotor, floris, 
                   windshear::Union{Matrix, WindShear})




=#


using LinearAlgebra
using StaticArrays
using Infiltrator

include("veer.jl")
include("structs_floris.jl")
include("runfloris.jl")
include("gaussian.jl")




buffers = FLORISBuffers(safehouse.nRP)
safehouse.set.enable_veer = true

#=
procedure for getting safehouse:

- uncomment @infiltrate in runFLORIS!
- include("examples/veer_main_mini.jl")
- [inside REPL] @exfiltrate to get safehouse
- remove @infiltrate again
- include("src/floris/testscripts.jl")    (this script)

now ready for runFLORIS dev
=#

#note, NOT safehouse.buffers
runFLORIS!(buffers, safehouse.set, safehouse.location_t, safehouse.states_wf, safehouse.states_t, safehouse.d_rotor, safehouse.floris, safehouse.windshear)   




















