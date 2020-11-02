#= 
アルゴリズム5 LLL簡約基底
=#
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
        end

        cond = lovasz_condition(δ, μ, B, k)
        if cond
            k += 1
        else
            b[k, :], b[k - 1, :] = b[k - 1, :], b[k, :]
            gso_update_lll!(μ, B, k)
            k = max(k - 1, 2)
        end
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