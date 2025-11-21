#funcs.jl


function teststruct(n_ops::Int)::teststruct
    return teststruct(
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Matrix{Float64}(undef, n_ops, n_ops),
        Matrix{Float64}(undef, n_ops, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops),
        Vector{Float64}(undef, n_ops)
    )
end





function setup_teststruct_buffers(n_ops::Int)::teststruct
    #...
end




function compute_wake_effects!(buf::teststruct, params::Dict{String,Any})
    #...
end