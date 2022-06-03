using Test, HopTB.Memoize


let
    ncalls = 0
    @memoize x function f(x, y::Int64=1; z::Int64=8)::Int64
        ncalls += 1
        return x + y + z
    end
    f(1, 2; z=3)
    @test ncalls == 1
    f(1, 2; z=3)
    @test ncalls == 1
    f(1, 3; z=4)
    @test ncalls == 2
    f(1, 2; z=3)
    @test ncalls == 2
    f(2, 3; z=5)
    @test ncalls == 3
    f(1, 2; z=3)
    @test ncalls == 4
    @test length(var"##_f_cache".data) == 1
end

let
    ncalls = 0
    @memoize function g(x, y::Int64=1; z::Int64=8)::Int64
        ncalls += 1
        return x + y + z
    end
    g(1, 2; z=3)
    @test ncalls == 1
    g(1, 2; z=3)
    @test ncalls == 1
    g(1, 3; z=4)
    @test ncalls == 2
    g(1, 2; z=3)
    @test ncalls == 2
    g(2, 3; z=5)
    @test ncalls == 3
    g(1, 2; z=3)
    @test ncalls == 3
    @test length(var"##_g_cache".data) == 3
end