#funcs.jl

#constructor for bufferstruct
function bufferstruct(nRP::Int)::bufferstruct
    return bufferstruct(
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Matrix{Float64}(undef, 2, 2),
        Matrix{Float64}(undef, 2, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),  
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP),
        Vector{Float64}(undef, nRP)
    )
end

function Params()::Params
    return Params(
        0.05,
    )
end





function setup_computation_buffers!(buf::bufferstruct, nRP::Int, nT::Int)
    buf.alpha .= 0
    buf.beta = Vector{Float64}(undef, nRP)
    buf.gamma = Vector{Float64}(undef, nRP)
    buf.a_star = Vector{Float64}(undef, nRP)
    buf.xi_0_hat = Vector{Float64}(undef, nRP)
    buf.rotmtx = Matrix{Float64}(undef, 2, 2)
    
end




function compute_wake_effects!(buf::bufferstruct, par::Params)
    #...
end










function runFUNCTIONS(buf::bufferstruct, par::Params, nRP::Int, nT::Int)
    setup_computation_buffers!(buf, nRP, nT)
    compute_wake_effects!(buf, par)
end