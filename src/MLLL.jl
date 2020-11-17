#= 
アルゴリズム9(P89) MLLL簡約基底(TODO)
縦長な行列 =#
function mlll_reduced!(b::Matrix{Int64}, δ::Float64)
    m, n = size(b)
    # g: GSO情報を計算している位置
    # k: LLL基底簡約を行っている位置
    # l: LLL基底簡約を行う最後の位置
    # z: 入れ替え先(後ろから入れ替える)

    z, g = m, 1
    # B: B′のノルムのメモ化
    B = zeros(Float64, m)
    μ = zeros(Float64, (m, m))
    b′ = zeros(Float64, (m, n))

    while g <= z
        if norm(b[g, :]) == 0
            if g < z
                b[g, :], b[z, :] = b[z, :], b[g, :]
            end
            z -= 1
        end

        b′[g, :] = b[g, :]
        # g行目のベクトルに対してグラムシュミットを実行する
        for j in 1:g - 1
            @assert B[j] != 0
            μ[g, j] = b[g, :] ⋅ b′[j, :] / B[j]
            b′[g, :] -= μ[g, j] * b′[j, :]
        end
        B[g] = norm(b′[g, :])^2
        μ[g, g] = 1
        # ==
        
        if g == 1 # 一回目はグラムシュミットする必要がない
            g = 2
            continue
        end
        
        l, k, stop = g, g, false

        _cnt = 0
        while k <= l && !stop
            b, μ = partial_size_reduce!(k, k - 1, b, μ)
            ν = μ[k, k - 1]
            @assert ν <= 1/2
            β = B[k] + ν^2 * B[k - 1]

            if β >= δ * B[k - 1] # Lovasz条件を満たしているか
                for j in k - 2:-1:1
                    b, μ = partial_size_reduce!(k, j, b, μ)
                end
                k = k + 1
            else
                if iszero(b[k, :])
                    if k < z
                        b[z, :], b[k, :] = b[k, :], b[z, :]
                    end
                    z -= 1
                    g = k
                    stop = true
                else # b_k != 0
                    b[k - 1, :], b[k, :] = b[k, :], b[k - 1, :]
                    for j in 1:k - 2
                        μ[k, j], μ[k - 1, j] = μ[k - 1, j], μ[k, j]
                    end

                    # テクい更新
                    if β != 0
                        # c_k-1 == b_k, c_k == b_k-1
                        # b_1 ... b_k-1, ck-1 は一次独立 -> GSOを更新できる
                        if B[k] == 0
                            # 一次従属性を取り除く
                            B[k - 1] = β
                            b′[k - 1, :] *= ν
                            μ[k, k - 1] = 1 / ν
                            for i in k + 1:l
                                μ[i, k - 1] /= ν
                            end
                        else # b_k == 0
                            # b_1 ... b_k-1 (k - 1こ)の一次従属性を考える
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
                        # println("B[$k] <-> B[$(k - 1)]")
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

        if !stop
            g += 1
        end
    end
end