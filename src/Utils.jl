#=
各基底ベクトル(行ベクトル)のノルム
=#
function hnorms(b::Matrix)::Array{Float64}
    [norm(b[i, :]) for i in 1:size(b, 1)]
end

function vol(b::Matrix)::Float64
    prod(hnorms(b))
end

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