@testset "DeepLLL" begin
    B = [
        9 2 7
        8 6 1
        3 2 6
    ]
    δ = 0.75
    B′ = copy(B)
    deep_lll_reduced!(B′, δ)
    B′_true = [
        3 -2 -5
        6 0 1
        2 6 0
    ]
    @test all(B′ == B′_true)
end