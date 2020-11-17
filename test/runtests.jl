using Test
using Base
using LinearAlgebra
using LatticeImpl
using DelimitedFiles

function Base.show(io::IO, m::Matrix{T}) where T
    println(io, "")
    writedlm(io, m)
    println(io, "")
end

function Base.show(io::IO, m::Matrix{Float64})
    writedlm(io, map(e -> round(e, digits=2), m))
    println(io, "")
end

function passed(f::Function, args...)::Bool
    try
        f(args...)
        true
    catch
        false
    end
end

@testset "All tests" begin
    include("./SizeReduce.jl")
    include("./LagrangeReduce.jl")
    include("./LLL.jl")
    include("./DeepLLL.jl")
    include("./MLLL.jl")
    include("./Babai.jl")
    include("./BDD.jl")
end