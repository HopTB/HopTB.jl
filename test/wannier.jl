include("zoo.jl")

let
    tm = getBN()
    twfs = Dict([0, 0, 0]=>reshape([1.0, 0.0], (2, 1)))
    wfs = Hop.Wannier.getwf(tm, twfs, [1], [10, 10, 1], [5, 5, 1])
    @test wfs[[0, 0, 0]][2] ≈ wfs[[-1, 0, 0]][2]
end

let
    tm = getBN()
    twfs = Dict([0, 0, 0]=>reshape([1.0, 0.0], (2, 1)))
    Rs = [0 0 0; -1 0 0]'
    wfs = Hop.Wannier.getwf(tm, twfs, [1], Rs, atol=1.0e-15)
    @test wfs[[0, 0, 0]][2] ≈ wfs[[-1, 0, 0]][2]
    nrmesh = [2, 2, 0]
    Rs = hcat([[Rx, Ry, Rz] for Rx in -nrmesh[1]:nrmesh[1] for Ry in -nrmesh[2]:nrmesh[2] for Rz in -nrmesh[3]:nrmesh[3]]...)
    wfs = Hop.Wannier.getwf(tm, twfs, [1], Rs)
    s = 0.0
    for (_, val) in wfs
        s += norm(val[:, 1])^2
    end
    @test s ≈ 1 atol=0.02
end

let
    # This is an example of constructing Wannier function of Liang Fu's C_4+T
    # TCI by interpolating the topological phase to atomic limit by breaking T.
    # The T symmetry is obtained by adding onsite energy to px±ipy orbitals on
    # different atoms.
    function getsTCI(t0, t1, t2, t3, t4, o1, o2, o3, o4)
        TCI = TBModel([1.0 0 0; 0 1 0; 0 0 1], copy([0 0 0; 0 0 0.2]'), [[1], [1]])
        ##################################################
        # original model
        ##################################################
        # blowup p_z orbitals
        sethopping!(TCI, [0, 0, 0], (1, 3), (1, 3), 10.0)
        sethopping!(TCI, [0, 0, 0], (2, 3), (2, 3), 10.0)
        # nearest neighbour
        addhopping!(TCI, [1, 0, 0], (1, 1), (1, 1), t0*2)
        addhopping!(TCI, [1, 0, 0], (2, 1), (2, 1), -t0*2)
        # next-nearest neighbour
        for i in 1:2, (n, m) in [(1, 1), (1, 2), (2, 1), (2, 2)]
            addhopping!(TCI, [1, 1, 0], (i, n), (i, m), (-1)^(i+1)*t1*2)
        end
        # interlayer orbital independent hopping
        for n in 1:2
            addhopping!(TCI, [0, 0, 0], (1, n), (2, n), t2)
            addhopping!(TCI, [0, 0, 1], (2, n), (1, n), t3)
            addhopping!(TCI, [1, 0, 0], (1, n), (2, n), t4*4)
        end
        ##################################################
        # symmetrization
        ##################################################
        C4 = Hop.Group.getrotation(2π/4, [0, 0, 1])
        nTCI = Hop.changebasis(TCI, Dict(1=>[-1/√2 im/√2 0; 0 0 1; 1/√2 im/√2 0]))
        set_is_canonical_ordered!(nTCI, true)
        sTCI = Hop.Group.symmetrize(nTCI, Hop.Group.generategroup([C4]))
        ##################################################
        # TRS breaking terms
        ##################################################
        addhopping!(sTCI, [0, 0, 0], (1, 1), (1, 1), o1)
        addhopping!(sTCI, [0, 0, 0], (1, 3), (1, 3), o2)
        addhopping!(sTCI, [0, 0, 0], (2, 1), (2, 1), o3)
        addhopping!(sTCI, [0, 0, 0], (2, 3), (2, 3), o4)
        ##################################################
        return sTCI
    end
    twfs = Dict([0, 0, 0]=>[0.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0])
    wfs = Hop.Wannier.interpolatewf(getsTCI, copy([
        0.0 0.0 0.0 0.0 0.0 1.0 -1.0 -1.0 1.0; # This is atomic limit
        0.0 0.0 0.0 0.0 0.5 1.0 -1.0 -1.0 1.0;
        0.0 0.0 0.0 2.0 0.5 1.0 -1.0 -1.0 1.0;
        0.0 0.0 2.5 2.0 0.5 1.0 -1.0 -1.0 1.0;
        0.0 0.25 2.5 2.0 0.5 1.0 -1.0 -1.0 1.0;
        1.0 0.25 2.5 2.0 0.5 1.0 -1.0 -1.0 1.0;
        1.0 0.25 2.5 2.0 0.5 0.0 0.0 0.0 0.0; # This is the parameter of the original model
    ]'), twfs, [1, 2], [10, 10, 10], [3, 3, 3], 3)
    @test abs(wfs[[0, 0, 0]][3, 1]) ≈ 0.6801019654070432
    @test wfs[[1, 0, 0]][1, 1] ≈ -wfs[[0, 1, 0]][1, 1]
end
