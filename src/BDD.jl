#= 
探索LWE問題の求解法(TODO)
探索LWE問題をBDDに帰着
BDDを解くdecoding attack
BDD: Bounded Distance Decoding

A: m * n ∈ Z_q^n
b: m * 1 ∈ Z_q 
秘密ベクトルとq_ary_baseの基底 =#
function decoding_attack(q::Int64, σ::Float64, A::Matrix{Int64}, b::Vector{Int64})::Vector{Int64}
    n = size(A', 1)
    m = size(A', 2)
    q_ary_base = vcat(A', q * Matrix{Int64}(I, m, m))
    @show size(q_ary_base)
    δ = 0.99
    q_ary_base = mlll_reduced!(q_ary_base, δ)
    @assert size(q_ary_base, 1) == size(q_ary_base, 2) == m
    v = babai_fast(q_ary_base, b)
    e = b - v
    A′ = A[n, :]
    b′ = b[n]
    e′ = e[n]
    s = (b′ - e′) ⋅ inv(A′')

    return s .% q
end