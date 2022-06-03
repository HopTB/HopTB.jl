using Hop, Test, LinearAlgebra
include("zoo.jl")

@register "Si.dat" nm->begin
    χ = Hop.Optics.getpermittivity(nm, 1, 1, collect(0:0.1:5), -5.87321, [10, 10, 10], ϵ=0.1)
    @test χ[40] ≈ -9.3976 + 12.2981im atol = 1.0e-3
end

# WS2
@register "WS2.dat" tm->begin
    @test Hop.Optics.get_shift_cond(tm, 2, 2, [3.2], -6.0, [5, 5, 1], ϵ=√0.1)[1] ≈ 7.7415 atol = 1.0e-4
    @test Hop.Optics.get_shift_cond(tm, 2, 1, [3.2], -6.0, [5, 5, 1], ϵ=√0.1)[1] ≈ -6.6951 atol = 1.0e-4
end

# BN model
let
    tm = Hop.Zoo.getBN()

    @test Hop.Optics.get_shift_cond(tm, 2, 2, [0.5, 1.0, 1.5], 0.0, [100, 100, 1], ϵ=√0.1) ≈
        [0.0391185, 1.3380897, 1.76826] atol = 1.0e-3
    
    @test Hop.Optics.get_shift_cond(tm, 2, 2, [0.5, 1.0, 1.5], 0.0, [100, 100, 1], ϵ=√0.1, batchsize=100) ≈
        [0.0391185, 1.3380897, 1.76826] atol = 1.0e-3
    
    @test Hop.Optics.get_shift_cond(tm, 2, 1, [0.5, 1.0, 1.5], 0.0, [100, 100, 1], ϵ=√0.1) ≈
        [-0.0391185, -1.3380897, -1.76826] atol = 1.0e-3
    
    @test Hop.Optics.get_shift_cond(tm, 2, 2, [0.5, 1.0, 1.5], 0.0, [10, 10, 1]; ϵ=√0.1, batchsize=13) ≈ 
        Hop.Optics.get_shift_cond(tm, 2, 2, [0.5, 1.0, 1.5], 0.0, [10, 10, 1]; ϵ=√0.1, batchsize=10)
    
    @test Hop.Optics.get_shift_cond_k(tm, 2, 2, 2, [3.6], 0.0, [1/3, 1/3, 0]; ϵ=0.1) ≈ [0.106791] atol=1.0e-3
    
    kpts = Hop.Utilities.constructlinekpts([0.0 1.0; 0.0 0.0; 0.0 0.0], tm.rlat, 5)[2]

    @test Hop.Optics.cltberryconnection(tm, 2, kpts)[1, 2, 3] ≈ 0.06667 atol = 1.0e-3

    @test Hop.Optics.cltshiftvector(tm, 2, 2, kpts)[1, 2, 2] ≈ 0.8819 atol = 1.0e-3
end

@register "WS2.w90" tm->begin
    @test Hop.Optics.getpermittivity(tm, 2, 2, [2.0,], -1.5, [10, 10, 1], ϵ=0.1)[1] ≈ 3.3628 + 0.7921im atol = 1.0e-3
end

let
    # stability test
    tm = get_Kane_Mele(include_rashba=true)
    @test real(Hop.Optics.get_shg(tm, 1, 2, 2, [0.5, 1.0, 1.5], 0.0, [5, 5, 1])[1]) ≈ 0.0014 atol=1.0e-4
    @test real(Hop.Optics.get_shg_B(tm, 1, 2, 2, [0.5, 1.0, 1.5], 0.0, 10.0, 3, [5, 5, 1])[1]) ≈ 0.06596 atol=1.0e-4
    tm = get_Kane_Mele(include_rashba=false)
    @test Hop.Optics.get_shg(tm, 1, 2, 2, [0.5, 1.0, 1.5], 0.0, [5, 5, 1])[1] ≈ 0.0 atol=1.0e-4
end

@register "GaAs.openmx" tm->begin
    @test Hop.Optics.get_shg(tm, 1, 2, 3, [1.5], -4.2, [5, 5, 5]; ϵ=0.1, scissor=0.69)[1] ≈ -255.88643+379.955523im atol=1.0e-3
end

@register "GaAs.openmx39" tm->begin
    @test Hop.Optics.get_shg(tm, 1, 2, 3, [1.5], -4.2, [5, 5, 5]; ϵ=0.1, scissor=0.69)[1] ≈ -255.84544+379.86693im atol=1.0e-3
end

let
    tm = get_Kane_Mele()
    @test Hop.Optics.get_Drude_weight(tm, 1, 1, -1.0, [10, 10, 1]; kBT=0.01) ≈ 761.9797953601502 + 0.0im
    @test Hop.Optics.get_Drude_weight(tm, 1, 2, -1.0, [10, 10, 1]; kBT=0.01) ≈ 0.0 atol=1.0e-8
end
