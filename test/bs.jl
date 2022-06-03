using Hop, Test, LinearAlgebra
include("zoo.jl")

@register "Si.dat" nm->begin
    kdist, egvals = Hop.BandStructure.getbs(nm, [0.0 0.5; 0.0 0.0; 0.0 0.0], 5)
    @test kdist[2] ≈ 0.25117 atol = 1.0e-3
    @test egvals[20, 1] ≈ 2.0243 atol = 1.0e-3
    @test Hop.BandStructure.clteig(nm, [0.0 0.5; 0.0 0.5; 0.0 0.5])[14, 2] ≈ -7.0692 atol = 1.0e-3
end

@register "WS2.scfout" nm->begin
    kdist, egvals = Hop.BandStructure.getbs(nm,
        [0 0 0;0.5 0 0;0.5 0 0;1 / 3 1 / 3 0;1 / 3 1 / 3 0;0 0 0]', 5)
    @test kdist[2] ≈ 0.28423 atol = 1.0e-3
    @test kdist[10] ≈ 1.79333 atol = 1.0e-3
    @test egvals[4, 3] ≈ -19.83008 atol = 1.0e-3
    @test egvals[10, 15] ≈ -7.65829 atol = 1.0e-3
end

@register "WS2-soc.scfout" nm->begin
    kdist, egvals = Hop.BandStructure.getbs(nm,
        [0 0 0;0.5 0 0;0.5 0 0;1 / 3 1 / 3 0;1 / 3 1 / 3 0;0 0 0]', 5)
    @test kdist[2] ≈ 0.28423 atol = 1.0e-3
    @test kdist[10] ≈ 1.79333 atol = 1.0e-3
    @test egvals[4, 3] ≈ -41.70918 atol = 1.0e-3
    @test egvals[10, 15] ≈ -19.16184 atol = 1.0e-3
end

@register "WS2.w90" nm->begin
    kdist, egvals = Hop.BandStructure.getbs(nm, copy([0 0 0; 0.5 0 0]'), 5)
    @test egvals[4, 3] ≈ -4.9114 atol = 1.0e-3
end

let
    # stability test
    tm = getBN()
    @test Hop.BandStructure.getjdos(tm, [1.5], 0.0, [10, 10, 1])[1] ≈ 0.01752601611101201
    @test Hop.BandStructure.clt_jdos(tm, [3.0], 0.0, [0.0 0.0 0.0; 0.6 0.6 0.0]')[2, 1] ≈ 0.4124189710369726
end


let
    tm = getBN()
    Es = Hop.BandStructure.clteig(tm, [5, 5, 1])
    @test geteig(tm, [0.2, 0.4, 0.0]).values ≈ Es[:, 2, 3, 1]
end

let
    tm = getBN()
    dos = Hop.BandStructure.getdos(tm, -5:0.01:5, [10, 10, 1], ϵ=0.01)
    @test sum(dos)*0.01 ≈ det(tm.rlat)*2/(2π)^3 rtol=1.0e-3
    @test dos[300] ≈ 0.25551954095619583 rtol=1.0e-3
end

@testset "Fermi surface extraction" begin
    tm = Hop.Zoo.getcube()
    fs = Hop.BandStructure.get_fermi_surfaces(tm, [50, 50, 50], [1]; fermi_energy=0.1)[1]
    @test all([isapprox(geteig(tm, fs.ks[:, i]).values[1], 0.1, atol=0.01) for i in 1:size(fs.ks, 2)])
end

@testset "EnergyMesh" begin
    tm = Hop.Zoo.getcube()
    em = Hop.BandStructure.get_energy_mesh(tm, [50, 50, 50], [1])
    nem = Hop.BandStructure.add_endboundary(em)
    @test nem.energies[:, end, end, end] ≈ em.energies[:, 1, 1, 1]
    @test norm(Hop.BandStructure.remove_endboundary(nem).energies - em.energies) < 1.0e-3
end

@testset "Density of states" begin
    tm = Hop.Zoo.getcube()
    em = Hop.BandStructure.get_energy_mesh(tm, [50, 50, 50], [1])
    dos = Hop.BandStructure.getdos(tm, em, -7:0.01:7; ϵ=0.1)
    @test sum(dos) * abs(det(tm.lat)) * 0.01 ≈ 1 atol=1.0e-3
end