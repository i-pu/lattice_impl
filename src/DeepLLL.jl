#= 
アルゴリズム7 DeepLLL簡約基底
=#
function deep_lll_reduced!(b::Matrix{Int64}, δ::Float64)
    @assert 1 / 4 < δ < 1
    n = size(b, 1)
    b′, μ = gram_schmidt(b)

    B = [norm(b′[j, :]) ^ 2 for j in 1:n]

    k = 2
    while k <= n
        for j in k - 1:-1:1
            partial_size_reduce!(k, j, b, μ)
        end

        C = norm(b[k, :]) ^ 2
        i = 1
        while i < k
            if C >= δ * B[i] 
                # |π_i(b_k)|^2を計算
                C = C - (μ[k, i]^2) * B[i]
                i += 1
            else
                v = b[k, :]            
                for j in k:-1:i+1
                    # 本のアルゴリズムには誤植がある
                    b[j, :] = b[j - 1, :]
                end
                b[i, :] = v

                # gso_update_deep_lll!(i, k, μ, B)
                # update μ, B
                b′, μ = gram_schmidt(b)
                B = [norm(b′[j, :]) ^ 2 for j in 1:n]

                k = max(i, 2) - 1
            end
        end
        k += 1
    end
end

#= 
DeepLLL GSO更新アルゴリズム
基底を交換(deep insertion)し効率的にμを再計算
=#
function gso_update_deep_lll!(i::Int64, k::Int64, μ::Matrix{Float64}, B::Array{Float64})
    @assert i < k
    n = size(μ, 1)
    P, D = copy(B), copy(B)
    for j in (k - 1):-1:i
        P[j] = μ[k, j] * B[j]
        # D[j] = ||π_j(b_k)||^2 = B_k + ∑_{k - 1}^{m = j} μ^2_k,m B_m
        D[j] = D[j + 1] + μ[k, j] * P[j]
    end

    S = zeros(size(B))

    # ths 2.4.3(1)
    for j in k:-1:(i + 1)
        T = μ[k, j - 1] / D[j]
        for l in n:-1:(k + 1)
            S[l] = S[l] + μ[l, j] * P[j]
            μ[l, j] = μ[l, j - 1] - T * S[l]
        end

        for l in k:-1:(j + 1)
            S[l] = S[l] + μ[l - 1, j] * P[j]
            μ[l, j] = μ[l - 1, j - 1] - T * S[l]
        end
    end

    # ths 2.4.3(2)
    T = 1 / D[i]
    for l in n:-1:(k + 1)
        μ[l, i] = T * (S[l] + μ[l, i] * P[i])
    end
    for l in k:-1:(i + 2)
        μ[l, i] = T * (S[l] + μ[l-1, i] * P[i])
    end
    μ[i + 1, i] = T * P[i]
    for j in 1:(i - 1)
        ϵ = μ[k, j]
        for l in k:-1:(i + 1)
            μ[l, j] = μ[l - 1, j]
        end
        μ[i, j] = ϵ
    end

    for j in k:-1:(i + 1)
        B[j] = D[j] * B[j - 1] / D[j - 1]
    end
    B[i] = D[i]
end