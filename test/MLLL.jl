@testset "MLLL" begin
    A = [
      1 5 21 3 14
      17 0 12 12 13
      12 21 15 6 6
      4 13 24 7 16
      20 9 22 27 8 
      19 8 19 3 1
      18 22 4 8 18
      10 11 19 18 21
      28 18 24 27 20
    ]
    B_true = [
      -2 0 -1 1 -4 1 -2 1 -1 0
      2 -4 0 0 1 1 0 -4 1 1
      3 -2 1 -1 1 0 0 3 3 -1
      -2 2 -3 2 0 3 0 0 -1 1
      2 4 -2 0 -3 0 0 0 0 -1
      -1 2 3 0 1 2 3 1 -2 0
      -1 -3 -1 -1 -2 -1 -1 3 0 4
      0 2 2 4 -2 -3 0 1 -2 3
      2 -1 -1 2 -1 1 1 5 -5 0
      -1 2 1 -2 3 3 -6 -3 0 1
    ]
    @assert size(B_true) == (10, 10)
    q, δ = 29, 0.99
    n = size(A', 1)
    m = size(A', 2)
    b = vcat(A', q * Matrix{Int64}(I, m, m))

    mlll_reduced!(b, δ)

    @show B
    @show B_true
    @test B == B_true
end