#= 
アルゴリズム5 LLL簡約基底 =#
using LinearAlgebra
include("./lagrange_reduced.jl")
include("./size_reduce.jl")


function lll_reduced!(b::Matrix{Int64}, δ::Float64)
    @assert 1 / 4 < δ < 1
    n = size(b, 1)
    b′, μ = gram_schmidt(b)

    B::Array{Float64} = []
    for i in 1:n
        push!(B, norm(b′[i, :]) ^ 2)
    end

    k = 2
    while k <= n
        for j in k - 1:-1:1
            partial_size_reduce!(k, j, b, μ)
            # println("after GSO ($k $j) で size_reduce")
            # @show b
            # @show μ
            # println("")
        end

        cond = lovasz_condition(δ, μ, B, k)
        # println("lovasz_condition $cond\n")
        if cond
            k += 1
        else
            b[k, :], b[k - 1, :] = b[k - 1, :], b[k, :]
            gso_update_lll!(μ, B, k)
            k = max(k - 1, 2)
        end

        # println("after $k")
        # @show b
        # @show μ
        # @show B
        # println("")
    end
end
    
#= 
GSO更新アルゴリズム
基底を交換し効率的にμを再計算
=#
function gso_update_lll!(μ::Matrix{Float64}, B::Array{Float64}, k::Int64)
    ν = μ[k, k - 1]
    n = size(μ, 1)
    b = B[k] + ν^2 * B[k - 1]
    # update 
    μ[k, k - 1] = ν * B[k - 1] / b
    B[k] = B[k] * B[k - 1] / b
    B[k - 1] = b
    
    for j in 1:(k - 2)
        μ[k, j], μ[k - 1, j] = μ[k - 1, j], μ[k, j]
    end

    for i in (k + 1):n
        t = μ[i, k]
        μ[i, k] = μ[i, k - 1] - ν * t
        μ[i, k - 1] = t + μ[k, k - 1] * μ[i, k]
    end
end

#= 
Lovasz条件
=#
function lovasz_condition(δ::Float64, μ::Matrix{Float64}, B::Array{Float64}, k::Int64)::Bool
    B[k] >= (δ - (μ[k, k - 1]^2)) * B[k - 1]
end

function lll_property1(b::Matrix{Int64}, δ::Float64)
    b′, μ = gram_schmidt(b)
    α = 4 / (4 * δ - 1)
    n = size(b, 1)
    for j in 1:n 
        for i in j:n
            @assert norm(b′[j, :])^2 <= (α^(i - j)) * (norm(b′[i, :])^2)
        end
    end
end

function lll_property2(b::Matrix{Int64}, δ::Float64)
    α = 4 / (4 * δ - 1)
    n = size(b, 1)
    @assert norm(b[1, :]) <= (α^((n - 1) / 4)) * (vol(b)^(1 / n))
end

function lll_property3(b::Matrix{Int64}, δ::Float64)
    # \lambda_i関数が作れないのでぱす
end

function lll_property4(b::Matrix{Int64}, δ::Float64)
    α = 4 / (4δ - 1)
    n = size(b, 1)
    @assert prod(norms(b)) <= (α ^ ((n - 1) / 4)) * vol(b)
end

function norms(b::Matrix{Int64})::Array{Float64}
    [norm(b[i, :]) for i in 1:size(b, 1)]
end

function vol(b::Matrix{Int64})::Float64
    prod(norms(b))
end

function test_lll()
    B = [
        9 2 7
        8 6 1
        3 2 6
    ]
    δ = 0.75
    lll_reduced!(B, δ)
    B_true = [
        -1 4 -6
        2 6 0
        3 -2 -5
    ]
    @assert B == B_true
    lll_property1(B, δ)
    lll_property2(B, δ)
    lll_property3(B, δ)
    lll_property4(B, δ)
end

test_lll()