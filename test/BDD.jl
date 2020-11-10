@testset "BDD" begin
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
  s_true = [
    7
    27
    14
    23
    26
  ]
  b = [
    28
    2
    24
    26
    11
    14
    7
    28
    27
    13
  ]
  e_true = [
    0
    0
    0
    0
    -1
    0
    0
    0
    0
    0
  ]
  v_true = [
    28
    2
    24
    16
    12
    14
    7
    28
    27
    13
  ]

  n, q, σ = 5, 29, 0.29
  s = decoding_attack(q, σ, A, b)
  @test s == s_true
end