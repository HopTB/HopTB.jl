using HopTB, Test
using LinearAlgebra
include("zoo.jl")

let
    tm_ref = get_Kane_Mele()
    tm = get_Kane_Mele()
    HopTB.Structure.permute_orbits!(tm, [2, 4, 1, 3])
    @test tm.hoppings[[0, 0, 0]][1, 2] ≈ tm_ref.hoppings[[0, 0, 0]][2, 4]
    @test tm.hoppings[[1, 0, 0]][3, 2] ≈ tm_ref.hoppings[[1, 0, 0]][1, 4]
    @test tm.positions[[0, 0, 0]][1][2, 2] ≈ tm_ref.positions[[0, 0, 0]][1][4, 4]
    @test !HopTB.has_full_information(tm)
end


let
    tm = get_Kane_Mele()
    slab = HopTB.Structure.makeslab(tm, 1, 2)
    @test slab.norbits == 8
    @test slab.orbital_types == [[0] for i in 1:4]
    @test slab.hoppings[[0, -1, 0]][1, 2] ≈ tm.hoppings[[0, -1, 0]][1, 2]
    @test slab.hoppings[[0, 0, 0]][3, 2] ≈ tm.hoppings[[-1, 0, 0]][1, 2]
    @test slab.hoppings[[0, 0, 0]][7, 6] ≈ tm.hoppings[[-1, 0, 0]][3, 4]
    @test slab.site_positions[1, 3] ≈ tm.positions[[0, 0, 0]][1][1, 1]+tm.lat[1, 1]
    @test slab.site_positions ≈ HopTB.get_orbital_position(slab)[:, 1:4]
    @test slab.positions[[0, 0, 0]][1][3, 3] ≈ tm.positions[[0, 0, 0]][1][1, 1]+tm.lat[1, 1]
end


let
    tm_ref = get_Kane_Mele()
    tm = get_Kane_Mele()
    HopTB.Structure.shift!(tm, [0.1, 0, 0])
    @test tm.site_positions ≈ HopTB.get_orbital_position(tm)[:, 1:2]
    @test tm.positions[[0, 0, 0]][1][1, 2] ≈ 0
end


let
    tm = get_Kane_Mele()
    @test HopTB.Structure.make_superlattice(tm, [1 0 0; -1 2 0; 0 0 1]')[1][:, 1] == [-1, 0, 0]
    sc = HopTB.Structure.make_supercell(tm, [1 0 0; -1 2 0; 0 0 1]')
    @test sc.hoppings[[0, 0, 0]][1, 3] ≈ tm.hoppings[[0, 1, 0]][2, 2]
    @test sc.hoppings[[0, 0, 0]][1, 4] ≈ conj(tm.hoppings[[-1, 0, 0]][1, 2])
    @test sc.hoppings[[-1, 0, 0]][1, 4] ≈ tm.hoppings[[0, 0, 0]][2, 1]
    @test sc.hoppings[[0, 0, 0]][5, 7] ≈ tm.hoppings[[0, 1, 0]][4, 4]
    @test HopTB.get_orbital_position(sc)[:, 1] ≈ HopTB.get_orbital_position(tm)[:, 2]-tm.lat[:, 1]
end


let
    tm = getBN()
    sctm = HopTB.Structure.make_supercell(tm, copy([1 0 0; -1 2 0; 0 0 1]'))
    slab = HopTB.Structure.makeslab(tm, 1, 2)
    @test geteig(sctm, [0.0,0.0,0.0]).values[1] ≈ -3.04138 atol = 1e-4
    @test geteig(slab, [0.0,0.0,0.0]).values[1] ≈ -2.609895 atol = 1e-4
end
