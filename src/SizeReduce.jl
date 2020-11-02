#= 
アルゴリズム4 サイズ簡約基底
=#
"""
GSO係数を1/2の範囲に収めるサイズ基底簡約
"""
function size_reduce!(b::Matrix{Int64})::Tuple{Matrix{Int64},Matrix{Float64}}
    n = size(b, 1)
    _, μ = gram_schmidt(b)
    c = b
    for i in 2:n
        for j in i - 1:-1:1 
            c, μ = partial_size_reduce!(i, j, c, μ)
        end
    end
    return c, μ
end

function partial_size_reduce!(i::Int64, j::Int64, b::Matrix{Int64}, μ::Matrix{Float64})::Tuple{Matrix{Int64},Matrix{Float64}}
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