#= 
アルゴリズム4 サイズ簡約基底
=#
using LinearAlgebra

"""
グラムシュミット
"""
function gram_schmidt(b::Matrix{Int64})::Tuple{Matrix{Float64},Matrix{Float64}}
    b′, μ = zeros(Float64, size(b)), Matrix{Float64}(I, size(b))

    for i in 1:size(b, 1)
        b′[i, :] = b[i, :]
        for j in 1:i - 1
            μ[i, j] = b[i, :]' * b′[j, :] / norm(b′[j, :]) ^ 2
            b′[i, :] = b′[i, :] - (μ[i, j] * b′[j, :])
        end
    end
    b′, μ
end

"""
GSO係数を1/2の範囲に収めるサイズ基底簡約
"""
function size_reduce!(b::Matrix{Int64})::Tuple{Matrix{Int64},Matrix{Float64}}
    n = size(b, 1)
    _, μ = gram_schmidt(b)
    c = b
    for i in 1:n
        for j in i - 1:-1:1 
            c, μ = size_reduce_helper!(i, j, c, μ)
        end
    end
    return c, μ
end

function size_reduce_helper!(i::Int64, j::Int64, b::Matrix{Int64}, μ::Matrix{Float64})::Tuple{Matrix{Int64},Matrix{Float64}}
    if abs(μ[i, j]) <= 1 / 2
        return b, μ
    end   
    q = round(μ[i, j])
    b[i, :] = b[i, :] - q * b[j, :]
    for l in 1:j
        μ[i, l] = μ[i, l] - q * μ[j, l]
    end
    return b, μ
end

function test()
    b = [
        5 -3 -7
        2 -7 -7
        3 -10 0
    ]
    b′, μ = gram_schmidt(b)
    μ_true = [
        1 0 0
        0.96385 1 0
        0.54216 1.3107 1
    ]
    @show μ, μ_true
    B, _ = size_reduce!(b)
    B_true = [
        5 -3 -7
        -3 -4 0
        1 -3 7
    ]
    @show B
    @show B_true
    b_sr′, μ_sr = gram_schmidt(B)
    @show μ_sr
end

test()