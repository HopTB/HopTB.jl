using HopTB, Test

let
    lat = [1 1 / 2 0; 0 sqrt(3) / 2 0; 0 0 1]
    tm = TBModel(2, lat)
    @test tm.rlat ≈ [2pi 0 0; -2pi / sqrt(3) 4pi / sqrt(3) 0; 0 0 2pi]
    @test_throws Exception sethopping!(tm, [0, 0, 0], [1, 1], 0.5im)
    sethopping!(tm, [0, 0, 0], 1, 1, 0.5)
    @test_throws Exception addhopping!(tm, [0, 0, 0], [1, 1], 0.5im)
    addhopping!(tm, [0, 0, 0], 1, 1, 0.5)
    addhopping!(tm, [0, 0, 1], 1, 2, 0.5im)
    @test HopTB.has_full_information(tm) == false
    @test tm.hoppings[[0, 0, 0]][1, 1] ≈ 1.0
    @test tm.hoppings[[0, 0, 1]][1, 2] ≈ 0.5im
    @test tm.hoppings[[0, 0, -1]][2, 1] ≈ -0.5im
    @test tm.overlaps[[0, 0, 0]][1, 1] ≈ 1.0
    setposition!(tm, [0, 0, 0], 1, 1, [0.5, 0.0, 0.0])
    set_orbital_types!(tm, [[0], [0]])
    @test tm.nsites == 2
    @test tm.site_norbits == [1, 1]
    @test tm.site_positions[:, 1] ≈ [0.5, 0.0, 0.0]
    @test tm.is_canonical_ordered == false
end


let
    lat = [1 1 / 2 0; 0 √3 / 2 0; 0 0 1]
    tm = TBModel(4, lat, isorthogonal = false)
    setposition!(tm, [0, 0, 0], 1, 1, [0.5, 0.0, 0.0])
    @test_throws Exception set_orbital_types!(tm, [[0], [0]])
    @test_throws Exception set_orbital_types!(tm, [[0], [0]], isspinful = true)
    @test ismissing(tm.nsites)
    setposition!(tm, [0, 0, 0], 3, 3, [0.5, 0.0, 0.0])
    set_orbital_types!(tm, [[0], [0]], isspinful = true)
    @test tm.nsites == 2
    @test tm.site_norbits == [2, 2]
    @test tm.site_positions[:, 1] ≈ [0.5, 0.0, 0.0]
    @test tm.is_canonical_ordered == false
    @test_throws Exception setposition!(tm, [0, 0, 0], 1, 1, [0, 0, 0])
    @test_throws Exception setoverlap!(tm, [0, 0, 0], 1, 1, 0.5im)
    setoverlap!(tm, [0, 0, 0], 1, 1, 0.5)
    @test tm.positions[[0, 0, 0]][1][1, 1] ≈ 0.25
    sethopping!(tm, [0, 1, 0], (1, 1), (2, 2), 1.0)
    @test tm.hoppings[[0, 1, 0]][1, 4] ≈ 1.0
    sethopping!(tm, [0, 1, 0], 1, 2, [0.5 0.0; 0.0 0.25])
    @test tm.hoppings[[0, 1, 0]][3, 4] ≈ 0.25
    @test tm.hoppings[[0, 1, 0]][1, 4] ≈ 0.0
    @test HopTB._to_orbital_index(tm, 1) == [1, 3]
end


let
    lat = [1 1 / 2 0; 0 √3 / 2 0; 0 0 1]
    site_positions = lat * ([1 / 3 1 / 3 0; 2 / 3 2 / 3 0]')
    tm = TBModel(lat, site_positions, [[0], [0]])
    @test HopTB._to_orbital_index(tm, (2, 1)) == 2
    @test HopTB._to_site_index(tm, 2) == (2, 1)
    sethopping!(tm, [0, 0, 0], 1, 1, 0.5)
    addhopping!(tm, [0, 0, 0], 1, 1, 0.5)
    addhopping!(tm, [0, 0, 1], 1, 2, 0.5im)
    addhopping!(tm, [1, 0, 0], (1, 1), (2, 1), 1.0)
    @test HopTB.has_full_information(tm) == true
    @test tm.hoppings[[0, 0, 0]][1, 1] ≈ 1.0
    @test tm.hoppings[[0, 0, 1]][1, 2] ≈ 0.5im
    @test tm.hoppings[[0, 0, -1]][2, 1] ≈ -0.5im
    @test tm.overlaps[[0, 0, 0]][1, 1] ≈ 1.0
    @test tm.hoppings[[1, 0, 0]][1, 2] ≈ 1.0
    @test tm.positions[[0, 0, 0]][1][1, 1] ≈ site_positions[1, 1]
end

let
    lat = [1 1 / 2 0; 0 √3 / 2 0; 0 0 1]
    site_positions = lat * ([1 / 3 1 / 3 0; 2 / 3 2 / 3 0]')
    tm = TBModel(lat, site_positions, [[0], [0]], isspinful = true)
    @test HopTB._to_orbital_index(tm, (2, 1)) == 2
    @test HopTB._to_orbital_index(tm, (2, 2)) == 4
    @test HopTB._to_site_index(tm, 2) == (2, 1)
    @test HopTB._to_site_index(tm, 3) == (1, 2)
    @test tm.norbits == 4
    @test tm.positions[[0, 0, 0]][1][1, 1] ≈ site_positions[1, 1]
    @test tm.positions[[0, 0, 0]][2][3, 3] ≈ site_positions[2, 1]
    addhopping!(tm, [1, 0, 0], (1, 1), (2, 2), 1.0)
    @test tm.hoppings[[1, 0, 0]][1, 4] ≈ 1.0
end

let
    lat = [1 1/2 0; 0 √3/2 0; 0 0 1]
    tm = TBModel(2, lat)
    addhopping!(tm, [0, 0, 0], [1 0.5; 0.5 1])
    addhopping!(tm, [1, 0, 0], [0.25 0.1im; 0.1im 0.25])
    @test tm.hoppings[[0, 0, 0]][1, 1] ≈ 1
    @test tm.hoppings[[0, 0, 0]][1, 2] ≈ 0.5
    @test tm.hoppings[[1, 0, 0]][1, 2] ≈ 0.1im
    @test tm.hoppings[[-1, 0, 0]][2, 1] ≈ -0.1im
    @test_throws Exception addhopping!(tm, [0, 0, 0], [1 0.25; 0.5 1])
    @test_throws Exception addhopping!(tm, [0, 0, 0], [1 0.25 0.0; 0.5 1 0.0])
end

let
    lat = [1 0 0; 1/2 √3/2 0; 0 0 1]'
    orbital_positions = lat*([1/3 1/3 0; 2/3 2/3 0]')
    tm = TBModel(lat, orbital_positions, isorthogonal=true)
    @test tm.overlaps[[0, 0, 0]][1, 1] ≈ 1.0
    @test tm.positions[[0, 0, 0]][2][1, 1] ≈ orbital_positions[2, 1]
    @test tm.positions[[0, 0, 0]][1][2, 2] ≈ orbital_positions[1, 2]
end
