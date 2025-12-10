#=
function runFLORIS!(buffers::FLORISBuffers, set::Settings, location_t, states_wf, states_t, d_rotor, floris, 
                   windshear::Union{Matrix, WindShear})




=#


using LinearAlgebra
using StaticArrays
using Infiltrator
using FLORIDyn




include("veer.jl")
include("structs_floris.jl")
include("runfloris.jl")
include("gaussian.jl")






runFLORIS!(safehouse.buffers::FLORISBuffers, safehouse.set, safehouse.location_t, safehouse.states_wf, safehouse.states_t, safehouse.d_rotor, safehouse.floris, safehouse.windshear)   



