using LinearAlgebra

include("zoo.jl")

@register "WS2.dat" nm->begin
    W = Hop.Topology.get_wilson_spectrum(nm, collect(1:53), [0 0 0; 1 0 0]', 100)
    @test W[2] ≈ -0.8678 atol = 1.0e-3
end

let
    haldane = getHaldane()
    @test Hop.Topology.get_wilson_spectrum(haldane, [1], [0 1 / 3 0; 1 1 / 3 0]', 1000)[1] / pi + 0.5497 ≈ 0.0 atol = 1.0e-3
    @test Hop.Topology.get_wilson_spectrum(haldane, [1], [0 1 / 6 0; 1 1 / 6 0]', 1000)[1] / pi + 0.8374 ≈ 0.0 atol = 1.0e-3
end

let
    tm = TBModel([1.0 0 0; 0 1 0; 0 0 1], reshape([0.5, 0.0, 0.0], (3, 1)), [[0]])
    @test mod(Hop.Topology.get_wilson_spectrum(tm, [1], [0.0 0.0 0.0; 1.0 0.0 0.0]', 10000)[1], 2π) ≈ π atol = 1.0e-3
end

let
    tm = getBN()
    function getv(k)
        return geteig(tm, k).vectors[:, 1]
    end
    kpts = zeros(3, 101)
    kpts[1, :] = 0.0:0.01:1.0
    v0 = geteig(tm, [0, 0, 0]).vectors[:, 1]
    gauge1 = Hop.Topology.get_smooth_gauge(tm, getv, v0, kpts, unwind=true)
    @test norm(gauge1[:, 1]-gauge1[:, 101]) < 1.0e-8
    for i in 1:20
        @test norm(gauge1[:, i]-gauge1[:, i+1]) < 0.02
    end
    gauge2 = Hop.Topology.get_smooth_gauge(tm, 1, [20, 20, 20])
end
