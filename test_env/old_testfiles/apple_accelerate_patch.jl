using AppleAccelerate
using LinearAlgebra
using Logging

# 1. FIX THE BLAS PATH
# We use the absolute path to the Accelerate framework wrapper
const ACCELERATE_PATH = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

try
    # Use the absolute path to ensure the linker finds it
    LinearAlgebra.BLAS.lbt_forward(ACCELERATE_PATH)
    @info "BLAS/LAPACK backend switched to Apple Accelerate successfully."
catch e
    @warn "Could not switch BLAS backend. Falling back to OpenBLAS."
    # Optional: print(e) to see the specific dlopen error
end

# 2. UPDATED MAPPINGS (Removed :pow which doesn't exist in Base)
const ACCEL_MAPPINGS = Dict(
    :exp   => AppleAccelerate.exp,
    :log   => AppleAccelerate.log,
    :sin   => AppleAccelerate.sin,
    :cos   => AppleAccelerate.cos,
    :sqrt  => AppleAccelerate.sqrt,
    :abs   => AppleAccelerate.abs,

    # Note: For x^y, use AppleAccelerate.pow directly if needed, 
    # but overloading Base.^ is technically complex.
)

for (orig, accel) in ACCEL_MAPPINGS
    @eval begin
        import Base: $orig
        # Apply to Vectors and Matrices for Float32 and Float64
        $orig(x::Union{Array{Float32}, Array{Float64}}) = $accel(x)
    end
end

@info "Vectorized math patches applied for Float32/64 Arrays."