#= 
アルゴリズム18, 19 Babai最近平面アルゴリズム
=#
function babai(b::Matrix{Int64}, w::Vector{Int64})::Vector{Int64}
    n = size(b, 1)
    b′, _ = gram_schmidt(b)
    w_pre = w
    y::Vector{Array{Int64, 1}} = []
    for i in n:-1:1
        l = w_pre ⋅ b′[i, :] / norm(b[i, :] )^ 2
        yi = round(l) * b[i, :]
        push!(y, yi)
        w_pre = w_pre - (l - round(l)) * b′[i, :] - yi
    end
    
    return sum(y)
end

function babai_fast(b::Matrix{Int64}, w::Vector{Int64})::Vector{Int64}
    n, m = size(b)
    @assert n == m "$n == $m"
    b′, _ = gram_schmidt(b)
    β = copy(w)
    for i in n:-1:1
        l = β ⋅ b′[i, :] / (norm(b′[i, :])^2)
        β -= round(l) * b[i, :]
    end
    return w - β
end