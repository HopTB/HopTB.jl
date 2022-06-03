using Hop, Test
include("zoo.jl")

let
    lat = [(√3) / 2 (√3) / 2 0; -1 / 2 1 / 2 0; 0 0 1]
    BN = TBModel(2, lat)
    sethopping!(BN, [0, 0, 0], 1, 1, -1 / 2)
    sethopping!(BN, [0, 0, 0], 2, 2, 1 / 2)
    sethopping!(BN, [0, 0, 0], 1, 2, -1)
    sethopping!(BN, [-1, 0, 0], 1, 2, -1)
    sethopping!(BN, [0, -1, 0], 1, 2, -1)
    setposition!(BN, [0, 0, 0], 1, 1, 1, (BN.lat * [1 / 3, 1 / 3, 0])[1])
    setposition!(BN, [0, 0, 0], 1, 1, 2, (BN.lat * [1 / 3, 1 / 3, 0])[2])
    setposition!(BN, [0, 0, 0], 1, 1, 3, (BN.lat * [1 / 3, 1 / 3, 0])[3])
    setposition!(BN, [0, 0, 0], 2, 2, 1, (BN.lat * [2 / 3, 2 / 3, 0])[1])
    setposition!(BN, [0, 0, 0], 2, 2, 2, (BN.lat * [2 / 3, 2 / 3, 0])[2])
    setposition!(BN, [0, 0, 0], 2, 2, 3, (BN.lat * [2 / 3, 2 / 3, 0])[3])

    @test Hop.Hall._getahc(BN, 1, 2, [1 / 3 2 / 3 0; 1 / 3 2 / 3 0]')[1] ≈ -3.0 atol = 1.0e-7
    @test Hop.Hall.getahc(BN, 1, 2, [4, 4, 1])[1] ≈ 0.0 atol = 1.0e-7
    @test Hop.Hall.collect_berry_curvature(BN, 1, 2, [1 / 3 2 / 3 0; 1 / 3 2 / 3 0]') ≈ [-1.5 -1.5; 1.5 1.5] atol = 1.0e-7
end

let
    tm = get_Kane_Mele()
    @test Hop.Hall.getshc(tm, 1, 2, 3, [10, 10, 1])[1, 1] ≈ -393.2616 atol=1.0e-3
end

@register "Fe.openmx" tm->begin
    ntm = Hop.change_energy_reference(tm, -5.441778786050317507)
    σs = Hop.Hall.getahc(ntm, 1, 2, [10, 10, 10])
    @test σs[1, 1] ≈ -794.815 atol=1.0e-3
end
