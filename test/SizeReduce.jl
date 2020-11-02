@testset "SizeReduce" begin
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
    @test isapprox(μ, μ_true; atol = 1e-3)
    B, _ = size_reduce!(b)
    B_true = [
        5 -3 -7
        -3 -4 0
        1 -3 7
    ]
    @test B == B_true
    b_sr′, μ_sr = gram_schmidt(B)
end