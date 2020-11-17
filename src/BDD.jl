#= 
探索LWE問題の求解法(TODO)
探索LWE問題をBDDに帰着
BDDを解くdecoding attack
BDD: Bounded Distance Decoding

A: m * n ∈ Z_q^n
b: m * 1 ∈ Z_q 
秘密ベクトルとq_ary_baseの基底 =#
using Mods
using LLLplus

function decoding_attack(q::Int64, σ::Float64, A::Matrix{Int64}, b::Array{Int64, 1})::Vector{Int64}
    m, n = size(A)
    @assert n <= m "$n <= $m"

    q_ary_base = [
        A' # m x n
        q * Matrix{Int64}(I, m, m) # m x m
    ]
    δ = 0.99999
    # q_ary_base, = seysen(q_ary_base)
    mlll_reduced!(q_ary_base, δ)
    @show q_ary_base
    @assert size(q_ary_base) == (n + m, m) "$(size(q_ary_base)) == $(n + m) x $m"
    @assert size(b, 1) == m "$(size(b, 1)) == $(m)"
    # v = babai_fast(q_ary_base[1:m, :], b)
    @show q_ary_base[1:m, :], Mod{q}
    v = babai(q_ary_base[1:m, :], b)
    v_true = [28, 2, 24, 16, 12, 14, 7, 28, 27, 13]
    @show norm(b - v), norm(b - v_true)
    @show v, v_true
    # v = [28, 2, 24, 16, 12, 14, 7, 28, 27, 13]
    e = b - v
    # [0, 0, 0, 0, -1, 0, 0, 0, 0, 0]
    A′ = Mod{q}.(A[1:n, :])
    b′ = b[1:n]
    e′ = e[1:n]
    @show size((b′ - e′)'), size(inv(A′')), inv(A′')
    @show b′, e
    @show A′' * inv(A′'), value.(inv(A′'))
    s = (b′ - e′)' * value.(inv(A′'))

    return s[1, :] .% q
end
