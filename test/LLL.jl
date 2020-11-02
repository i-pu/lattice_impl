function lll_property1(b::Matrix{Int64}, δ::Float64)
    b′, μ = gram_schmidt(b)
    α = 4 / (4 * δ - 1)
    n = size(b, 1)
    for j in 1:n 
        for i in j:n
            @assert norm(b′[j, :])^2 <= (α^(i - j)) * (norm(b′[i, :])^2)
        end
    end
end

function lll_property2(b::Matrix{Int64}, δ::Float64)
    α = 4 / (4 * δ - 1)
    n = size(b, 1)
    @assert norm(b[1, :]) <= (α^((n - 1) / 4)) * (vol(b)^(1 / n))
end

function lll_property3(b::Matrix{Int64}, δ::Float64)
    # \lambda_i関数が作れないのでぱす
end

function lll_property4(b::Matrix{Int64}, δ::Float64)
    α = 4 / (4δ - 1)
    n = size(b, 1)
    @assert prod(hnorms(b)) <= (α ^ ((n - 1) / 4)) * vol(b)
end

@testset "LLL" begin
    B = [
        9 2 7
        8 6 1
        3 2 6
    ]
    δ = 0.75
    lll_reduced!(B, δ)
    B_true = [
        -1 4 -6
        2 6 0
        3 -2 -5
    ]
    @test B == B_true
    @test passed(lll_property1, B, δ)
    @test passed(lll_property2, B, δ)
    @test passed(lll_property3, B, δ)
    @test passed(lll_property4, B, δ)
end