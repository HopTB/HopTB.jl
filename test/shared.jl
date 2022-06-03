using Hop, Test

let
    tm = Hop.Interface.createmodelopenmx(joinpath(dirname(@__FILE__), "data", "GaAs.openmx39"))
    sm = SharedTBModel(tm)
    tm_result = Hop.Optics.get_shg(tm, 1, 2, 3, [1.5], -4.2, [5, 5, 5]; ϵ=0.1, scissor=0.69)[1]
    sm_result = Hop.Optics.get_shg(sm, 1, 2, 3, [1.5], -4.2, [5, 5, 5]; ϵ=0.1, scissor=0.69)[1]
    @test tm_result ≈ sm_result atol=1.0e-4
end