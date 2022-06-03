
let
    @test Hop.Utilities.integrate(x->sin(x[1]+x[2]), [0.0, 0.0],
        [1.0, 1.0])[1] ≈ 2*sin(1)-sin(2)
    @test Hop.Utilities.integrate(x->sin(x[1]+x[2]), [0.0, 0.0, 0.0],
        [1.0, 1.0, 2.0], constant_components=[3])[1] ≈ 4*sin(1)-2*sin(2)
    @test Hop.Utilities.pintegrate(x->sin(x[1]+x[2]), [0.0, 0.0, 0.0],
        [1.0, 1.0, 2.0], constant_components=[3])[1] ≈ 4*sin(1)-2*sin(2)
    @test Hop.Utilities.pintegrate(x->sin(x[1]+x[2]), [0.0, 0.0],
        [1.0, 1.0])[1] ≈ 2*sin(1)-sin(2)
    @test Hop.Utilities.pintegrate(x->sin(x[1]+x[2]), [0.0, 0.0, 0.0],
        [1.0, 1.0, 2.0], constant_components=[3], ndiv=2)[1] ≈ 4*sin(1)-2*sin(2)
end

let
    @test Hop.Utilities.distance_on_circle(0.9, 0.1) ≈ 0.2
end


let
    lat = [
      0.     3.245  3.245;
      3.245  0.     3.245;
      3.245  3.245  0.;
    ]
    rlat = 2π*inv(lat)'
    kdist, kpts = Hop.Utilities.constructlinekpts([0.375 0.0 0.0 0.5; 0.375 0.0 0.0 0.5; 0.75 0.0 0.0 0.5], rlat, 10)
    @test kpts[1, 1] ≈ 0.375
end
