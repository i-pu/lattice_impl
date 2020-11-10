#= 
アルゴリズム9(P89) MLLL簡約基底(TODO) =#
function mlll_reduced!(b::Matrix{Int64}, δ::Float64)
    h = size(b, 1)
    # g: GSO情報を計算している位置
    # k: LLL基底簡約を行っている位置
    # l: LLL基底簡約を行う最後の位置
    z, g = h, 1
    while g <= z
        if iszero(b[g, :])
            if g < z
                b[g, :], b[z, :] = b[z, :], b[g, :]
            end
            z -= 1
        end

        # == GSOベクトルの計算
        b′, μ = similar(b, Float64), zeros((h, h))
        b′[g, :] = b[g, :]
        B = zeros(Float64, g)
        B[g] = norm(b[g, :])^2
        for j in 1:g-1
            # ここが0になってしまう
            @assert B[j] != 0
            μ[g, j] = b[g, :] ⋅ b′[j, :] / B[j]
            b′[g, :] -= μ[g, j] * b′[j, :]
        end
        B[g] = norm(b′[g, :]) ^ 2
        μ[g, g] = 1
        # ==
        
        if g == 1
            g = 2
        else
            l, k, start_again = g, g, false
            while k <= l && !start_again
                b, μ = partial_size_reduce!(k, k - 1, b, μ)
                ν = μ[k, k - 1]
                β = B[k] + ν^2 * B[k - 1]
                @show β >= δ * B[k - 1]
                if β >= δ * B[k - 1]
                    for j in k - 2:-1:1
                        b, μ = partial_size_reduce!(k, j, b, μ)
                    end
                    k = k + 1
                    if iszero(b[k, :])
                        if k < z
                            b[z, :], b[k, :] = b[k, :], b[z, :]
                        end
                        z -= 1
                        g = k
                        start_again = true
                        print("startagain")
                    else # b_k == 0
                        b[k - 1, :], b[k, :] = b[k, :], b[k - 1, :]
                        for j in 1:k - 2
                            μ[k, j], μ[k - 1, j] = μ[k - 1, j], μ[k, j]
                        end

                        if β != 0
                            if B[k] == 0
                                B[k - 1] = β
                                b′[k - 1, :] *= ν
                                μ[k, k - 1] = 1 / ν
                                for i in k + 1:l
                                    μ[i, k - 1] /= ν
                                end
                            else # b_k == 0
                                t = B[k - 1] / β
                                μ[k, k - 1] = ν * t
                                w = b′[k - 1, :]
                                b′[k - 1, :] = b′[k, :] + ν * w
                                B[k - 1] = β
                                if k <= l
                                    b′[k, :] = - μ[k, k - 1] * b′[k, :] + (B[k] / β) * w
                                    B[k] *= t
                                end

                                for i in k + 1:l
                                    t = μ[i, k]
                                    μ[i, k] = μ[i, k - 1] - ν * t
                                    μ[i, k - 1] = t + μ[k, k - 1] * μ[i, k]
                                end
                            end
                        else # β != 0
                            B[k], B[k - 1] = B[k - 1], B[k]
                            b′[k, :], b′[k - 1, :] = b′[k - 1, :], b′[k, :]
                            for i in k + 1:l
                                μ[i, k], μ[i, k - 1] = μ[i, k - 1], μ[i, k]
                            end
                        end

                        k = max(k - 1, 2)
                    end
                end
            end

            if !start_again
                g += 1
            end
        end
    end
end