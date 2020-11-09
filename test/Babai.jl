@testset "Babai" begin
    B = [
        1 2 3
        3 0 -3
        3 -7 3
    ]
    w = [
        10
        6
        5
    ]
    v = babai(B, w)
    v_true = [
        10
        8
        6
    ]
    @test v == v_true

    B = [
        1 2 3
        3 0 -3
        3 -7 3
    ]
    w = [
        10
        6
        5
    ]
    v = babai_fast(B, w)
    v_true = [
        10
        8
        6
    ]
    @test v == v_true
end

@testset "Babai Profile" begin
    B = [
        1 2 3
        3 0 -3
        3 -7 3
    ]
    w = [
        10
        6
        5
    ]
    @time for i in 0:100000; babai(B, w); end
    @time for i in 0:100000; babai_fast(B, w); end
end