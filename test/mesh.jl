using HopTB, Test


@testset "Meshes" begin
    # endboundary = false
    mesh = HopTB.Meshes.UniformMesh([3, 4, 5]; startpoint=[0.1, 0.2, 0.3], regionsize=[1.0, 2.0, 3.0])
    @test mesh[29] ≈ [1/3, 1/4, 2/5] .* [1.0, 2.0, 3.0] + [0.1, 0.2, 0.3]
    submesh = mesh[20:29]
    @test length(submesh) == 10
    @test submesh[10] ≈ [1/3, 1/4, 2/5] .* [1.0, 2.0, 3.0] + [0.1, 0.2, 0.3]
    @test HopTB.Meshes.find_point_in_mesh(mesh[29], mesh) == [2, 2, 3]
    # endboundary = true
    mesh = HopTB.Meshes.UniformMesh([3, 4, 5];
        startpoint=[0.1, 0.2, 0.3], regionsize=[1.0, 2.0, 3.0], endboundary=true)
    @test mesh[length(mesh)] ≈ [0.1, 0.2, 0.3] + [1.0, 2.0, 3.0]
    @test mesh[29] ≈ [1/2, 1/3, 2/4] .* [1.0, 2.0, 3.0] + [0.1, 0.2, 0.3]
    submesh = mesh[20:29]
    @test length(submesh) == 10
    @test submesh[10] ≈ [1/2, 1/3, 2/4] .* [1.0, 2.0, 3.0] + [0.1, 0.2, 0.3]
    @test HopTB.Meshes.find_point_in_mesh(mesh[29], mesh) == [2, 2, 3]
end