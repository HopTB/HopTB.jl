using Hop, Test, LinearAlgebra
include("zoo.jl")


@register "Si.dat" (nm->begin
    egvals = geteig(nm, [0.5, 0.5, 0.5]).values
    @test egvals[18] ≈ 4.5942 atol = 1.0e-3
end)

let
    tm = get_Kane_Mele()
    @test real.(diag(Hop.getvelocity(tm, 1, inv(tm.rlat)*[0.1, 0.2, 0.0]))) ≈
        (geteig(tm, inv(tm.rlat)*[0.1001, 0.2, 0.0]).values-geteig(tm, inv(tm.rlat)*[0.1, 0.2, 0.0]).values)/0.0001 atol=1.0e-3
    egvals = geteig(tm, [0.1, 0.2, 0.0]).values
    @test Hop.getvelocity(tm, 1, [0.1, 0.2, 0.0])[1, 2] ≈ im*Hop.getA(tm, 1, [0.1, 0.2, 0.0])[1, 2]*(egvals[1]-egvals[2])
end

@testset "Berry curvature" begin
    tm = Hop.Zoo.getBN()
    @test get_berry_curvature(tm, 1, 2, [1/3, 2/3, 0]) ≈ [-1.5, 1.5] atol = 1.0e-7
    @test get_berry_curvature(tm, 1, 2, [2/3, 1/3, 0]) ≈ [1.5, -1.5] atol = 1.0e-7
end
