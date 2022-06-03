using LinearAlgebra, Test, Hop

let
    C = Hop.Group.getrotation(pi / 4, [0, 0, 1], isdouble = true)
    @test Hop.Group.decompose_to_primitive(C).θ ≈ pi / 4
end

let
    T = Hop.Group.getTRS()
    @test T == T * T * T
    T = Hop.Group.getTRS(isdouble = true)
    @test !(T == T * T * T)
    C = Hop.Group.getrotation(2pi, [0, 0, 1], isdouble = true)
    @test T * T == C
end

let
    C = Hop.Group.getrotation(pi / 2, [0, 0, 1.0])
    @test Hop.Group.decompose_to_primitive(Hop.Group.inverse(C))[:θ] ≈ -pi / 2
    @test Hop.Group.decompose_to_primitive(Hop.Group.inverse(C))[:Ē] == false
    C = Hop.Group.getrotation(pi / 2 + 2pi, [0, 0, 1.0])
    @test Hop.Group.decompose_to_primitive(Hop.Group.inverse(C))[:θ] ≈ -pi / 2
    @test Hop.Group.decompose_to_primitive(Hop.Group.inverse(C))[:Ē] == true
    T = Hop.Group.getTRS()
    @test Hop.Group.decompose_to_primitive(Hop.Group.inverse(C * T))[:θ] ≈ -pi / 2
end

let
    Hop.Group.getJ(1)[2][1, 2] ≈ -im / √2
    T = Hop.Group.getTRS(isdouble = true)
    @test Hop.Group.get_orb_rep(T, 0) ≈ -im * [0.0 -im; im 0.0]
    @test Hop.Group.get_orb_rep(T, [0, 0])[1:2, 3:4] ≈ [-1 0; 0 -1]
    T = Hop.Group.getTRS(isdouble = false)
    @test Hop.Group.get_orb_rep(T, 1) ≈ [0 0 -1; 0 1 0; -1 0 0]
    C = Hop.Group.getrotation(pi, [0, 0, 1])
    @test C.rotation_matrix ≈ [-1.0 0 0; 0 -1 0; 0 0 1]
end

let
    lat = [(√3) / 2 (√3) / 2 0; -1 / 2 1 / 2 0; 0 0 1]
    site_positions = lat * ([1 / 3 1 / 3 0; 2 / 3 2 / 3 0]')
    orbital_types = [[0], [0]]
    tm = TBModel(lat, site_positions, orbital_types)
    sethopping!(tm, [0, 0, 0], 1, 2, 1.0)
    C3 = Hop.Group.getrotation(2pi / 3, [0, 0, 1], isdouble = false)
    @test Hop.Group._get_transformed_site(C3, tm, [0, 0, 0], 1)[1] ≈ [-1, 0, 0]
    @test Hop.Group._get_transformed_site(C3, tm, [0, 0, 0], 1)[2] ≈ 1
    set_is_canonical_ordered!(tm, true)
    stm = Hop.Group.symmetrize(tm, Hop.Group.generategroup([C3]))
    @test geteig(stm, [0.1, 0.0, 0.0]).values ≈ geteig(stm, [0.0, 0.1, 0.0]).values
end

let
    C3 = Hop.Group.getrotation(2pi / 3, [0, 0, 1], isdouble = true)
    @test length(Hop.Group.generategroup([C3])) == 6
    C3 = Hop.Group.getrotation(2pi / 3, [0, 0, 1], isdouble = false)
    @test length(Hop.Group.generategroup([C3])) == 3
    C3 = Hop.Group.getrotation(2pi / 3, [0, 0, 1], isdouble = true)
    My = Hop.Group.getmirror([0, 1, 0], isdouble = true)
    @test length(Hop.Group.generategroup([C3, My])) == 12
end

let
    lat = [(√3) / 2 (√3) / 2 0; -1 / 2 1 / 2 0; 0 0 1]
    site_positions = lat * ([1 / 3 1 / 3 0; 2 / 3 2 / 3 0]')
    orbital_types = [[1], [1]]
    tm = TBModel(lat, site_positions, orbital_types, isspinful = true)
    lat = [(√3) / 2 (√3) / 2 0; -1 / 2 1 / 2 0; 0 0 1]
    for n in [2, 5, 8, 11] # suppress pz orbitals
        sethopping!(tm, [0, 0, 0], n, n, -10.0)
    end
    sethopping!(tm, [0, 0, 0], 1, 1, 1.0)
    sethopping!(tm, [0, 0, 0], 3, 3, -0.5)
    sethopping!(tm, [0, 0, 0], 1, 4, -1.0 + 0.1im)
    sethopping!(tm, [0, 0, 0], 3, 6, -0.5 + 0.5im)
    sethopping!(tm, [0, 0, 0], 1, 10, 0.1)
    sethopping!(tm, [0, 0, 0], 1, 12, 0.3)

    C3 = Hop.Group.getrotation(2pi / 3, [0, 0, 1], isdouble = true)
    My = Hop.Group.getmirror([0, 1, 0], isdouble = true)
    set_is_canonical_ordered!(tm, true)
    stm = Hop.Group.symmetrize(tm, Hop.Group.generategroup([C3, My]))
    @test gethopping(stm, [0, 0, 0], (1, 1), (2, 4)) ≈ -gethopping(stm, [0, 0, 0], (1, 6), (2, 3))
    @test angle(gethopping(stm, [-1, 0, 0], (1, 1), (2, 4)) / gethopping(stm, [0, 0, 0], (1, 1), (2, 4))) ≈ -pi * 2 / 3
    @test gethopping(stm, [0, 0, 0], (1, 1), (2, 6)) ≈ gethopping(stm, [-1, 0, 0], (1, 1), (2, 6))
end

let
    lat = [(√3) / 2 (√3) / 2 0; -1 / 2 1 / 2 0; 0 0 1]
    site_positions = lat * ([1 / 3 1 / 3 0; 2 / 3 2 / 3 0]')
    orbital_types = [[0], [0]]
    tm = TBModel(lat, site_positions, orbital_types)
    sethopping!(tm, [0, 0, 0], 1, 2, 1.0)
    sethopping!(tm, [0, 0, 0], 1, 1, 1.0)
    C3 = Hop.Group.getrotation(2pi / 3, [0, 0, 1], isdouble = false)
    set_is_canonical_ordered!(tm, true)
    stm = Hop.Group.symmetrize(tm, Hop.Group.generategroup([C3]))
    @test stm.positions[[0, 0, 0]][1][1, 1] ≈ tm.positions[[0, 0, 0]][1][1, 1]
    @test Hop.Group.get_bloch_rep(C3, stm, [1 / 2, 1 / 2, 0]) ≈ [-1 0; 0 1]
end

let
    lat = Hop.Utilities.eye(3)
    site_positions = copy([0.2 0.0 0.0; 0.2 0.0 0.1; -0.2 0.0 0.0; -0.2 0.0 -0.1]')
    orbital_types = [[1], [0], [1], [0]]
    tm = TBModel(lat, site_positions, orbital_types, isspinful = false, is_canonical_ordered = true)
    sethopping!(tm, [0, 0, 0], (1, 2), (2, 1), 1.0)
    i = Hop.Group.getinversion(isdouble = false)
    stm = Hop.Group.symmetrize(tm, Hop.Group.generategroup([i]))
    @test gethopping(stm, [0, 0, 0], (1, 2), (2, 1)) ≈ 0.5
    @test gethopping(stm, [0, 0, 0], (3, 2), (4, 1)) ≈ -0.5
end

let
    lat = Hop.Utilities.eye(3)
    site_positions = copy([0.2 0.0 0.0; 0.2 0.0 0.1; -0.2 0.0 0.0; -0.2 0.0 -0.1]')
    orbital_types = [[1], [0], [1], [0]]
    tm = TBModel(lat, site_positions, orbital_types, isspinful = false, is_canonical_ordered = true)
    setposition!(tm, [0, 0, 0], 1, 4, 3, 1.0)
    i = Hop.Group.getinversion(isdouble = false)
    stm = Hop.Group.symmetrize(tm, Hop.Group.generategroup([i]))
    @test stm.positions[[0, 0, 0]][3][1, 4] ≈ stm.positions[[0, 0, 0]][3][5, 8]
    lat = Hop.Utilities.eye(3)
    site_positions = copy([0.2 0.0 0.0; 0.2 0.0 0.1; -0.2 0.0 0.0; -0.2 0.0 -0.1]')
    orbital_types = [[0], [0], [0], [0]]
    tm = TBModel(lat, site_positions, orbital_types, isspinful = false, is_canonical_ordered = true)
    setposition!(tm, [0, 0, 0], 1, 2, 3, 1.0)
    i = Hop.Group.getinversion(isdouble = false)
    stm = Hop.Group.symmetrize(tm, Hop.Group.generategroup([i]))
    @test stm.positions[[0, 0, 0]][3][1, 2] ≈ -stm.positions[[0, 0, 0]][3][3, 4]
end

let
    lat = convert(Array{Float64}, I(3))
    site_positions = lat*([0 -1/5 0; 0 -2/5 0; 1/2 1/5 0; 1/2 2/5 0]')
    orbital_types = [[0], [0], [0], [0]]
    tm = TBModel(lat, site_positions, orbital_types, isorthogonal=false)
    tm.is_canonical_ordered = true
    sethopping!(tm, [0, 0, 0], 1, 1, 1.0)
    setoverlap!(tm, [0, 0, 0], 1, 2, 0.1)
    setposition!(tm, [0, 0, 0], 1, 2, 2, 0.1)
    s = Hop.Group.gettranslation([1/2, 0, 0])*Hop.Group.getrotation(π, [1, 0, 0])
    Hop.Group.generategroup([s])
    ntm = Hop.Group.symmetrize(tm, Hop.Group.generategroup([s]))
    @test ntm.hoppings[[0, 0, 0]][1, 1] ≈ ntm.hoppings[[0, 0, 0]][3, 3]
    @test ntm.overlaps[[0, 0, 0]][1, 2] ≈ ntm.overlaps[[0, 0, 0]][3, 4]
    @test ntm.positions[[0, 0, 0]][2][1, 2] ≈ -ntm.positions[[0, 0, 0]][2][3, 4]
end

@register "SnH.openmx" SnH->begin
    orbital_types = [[0, 0, 1, 1, 2, 2, 2], [0, 0, 1, 1, 2, 2, 2], [0, 0, 1], [0, 0, 1]]
    set_orbital_types!(SnH, orbital_types, isspinful = true)
    kdist, egvals = Hop.BandStructure.getbs(SnH, [0.5 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 1 / 3 1 / 3 0.0]', 100)
    C3 = Hop.Group.getrotation(2pi / 3, [0, 0, 1], isdouble = true)
    Mx = Hop.Group.getmirror([1, 0, 0], isdouble = true)
    T = Hop.Group.getTRS(isdouble = true)
    i = Hop.Group.getinversion(isdouble = true)
    nSnH = changebasis(SnH, Hop.Group.Us_openmx)
    set_is_canonical_ordered!(nSnH, true)
    sSnH = Hop.Group.symmetrize(nSnH, Hop.Group.generategroup([C3, Mx, T, i]))
    Hop.prune!(sSnH)
    skdist, segvals = Hop.BandStructure.getbs(sSnH, [0.5 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 1 / 3 1 / 3 0.0]', 100)
    @test maximum(abs.(segvals - egvals)) < 1.0e-2
    C3rep = Hop.Group.get_bloch_rep(C3, sSnH, [0, 0, 0])
    Mxrep = Hop.Group.get_bloch_rep(Mx, sSnH, [0, 0, 0])
    v = geteig(sSnH, [0, 0, 0]).vectors[:, 29:30]
    @test eigvals(v' * getS(sSnH, [0, 0, 0]) * C3rep * v) ≈ [-1.0, -1.0]
    @test sort(eigvals(v' * getS(sSnH, [0, 0, 0]) * Mxrep * v), by = imag) ≈ [-im, im]
end

@register "WS2.openmx" WS2->begin
    orbital_types = [[0, 0, 0, 1, 1, 1, 2, 2, 3], [0, 0, 0, 1, 1, 1, 2, 2, 3], [0, 0, 0, 1, 1, 2, 2, 3]]
    set_orbital_types!(WS2, orbital_types)
    WS2 = convert(TBModel{ComplexF64}, WS2)
    WS2 = changebasis(WS2, Hop.Group.Us_openmx)
    set_is_canonical_ordered!(WS2, true)
    C3 = Hop.Group.getrotation(2π/3, [0, 0, 1])
    My = Hop.Group.getmirror([1, 0, 0])
    T = Hop.Group.getTRS()
    sWS2 = Hop.Group.symmetrize(WS2, Hop.Group.generategroup([C3, My, T]))
    for R in keys(WS2.positions)
        for α in 1:3
            @test maximum(abs.(WS2.positions[R][α]-sWS2.positions[R][α])) < 1.0e-4
        end
    end
end
