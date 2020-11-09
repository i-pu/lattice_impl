using Test
using Base
using LinearAlgebra
using LatticeImpl
using DelimitedFiles

function Base.show(io::IO, m::Matrix{T}) where T
    print(io)
    writedlm(stdout, m)
end

function Base.show(io::IO, m::Matrix{Float64})
    print(io)
    writedlm(stdout, map(e -> round(e, digits=2), m))
end

function passed(f::Function, args...)::Bool
    try
        f(args...)
        true
    catch
        false
    end
end

include("./SizeReduce.jl")
include("./LagrangeReduce.jl")
include("./LLL.jl")
include("./DeepLLL.jl")
include("./MLLL.jl")
include("./Babai.jl")