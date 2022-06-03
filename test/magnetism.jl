using Hop, Test, LinearAlgebra

include("zoo.jl")

let
    tm = get_Kane_Mele()
    @test Hop.Magnetism.get_orbital_moment(tm, 3, [0.1, 0.2, 0.0])[2, 2] ≈ 0.0074738780301569566
    @test Hop.Magnetism.get_field_modified_Es(tm, 3, 1.0, [0.1, 0.2, 0.0])[2] ≈ -2.5671018333590077
end
