@testset "LagrangeReduce" begin
  b1, b2 = [-7, -4, -10], [9, 5, 12]
  b1, b2 = ragrange_reduced(b1, b2)
  @test b1 == [1, 0, -2] && b2 == [2, 1, 2]

  b1, b2 = [230, -651, 609, -366], [301, -852, 797, -479]
  b1, b2 = ragrange_reduced(b1, b2)
  @test b1 == [-1, -3, -2, -1] && b2 == [2, -3, 5, -2]
end