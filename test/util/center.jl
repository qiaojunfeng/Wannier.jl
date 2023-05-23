using Wannier:Vec3
@testset "find_nearest_atom" begin
    wout = read_wout(joinpath(FIXTURE_PATH, "valence/band/silicon.wout"))

    points = map(c -> inv(wout.lattice) * c, wout.centers)  # to fractional
    distances, indexes, translations = Wannier.find_nearest_atom(
        points, wout.lattice, wout.atom_positions
    )

    ref_distances = [1.15401, 1.16433, 1.15363, 1.14736]
    ref_indexes = [2, 1, 1, 1]
    ref_translations = [
        0 0 0 1
        0 0 -1 -1
        0 0 0 0
    ]
    ref_translations = [Vec3(ref_translations[:, i]) for i = 1:size(ref_translations,2)]

    @test all(isapprox.(distances, ref_distances; atol=2e-5))
    @test indexes == ref_indexes
    @test translations == ref_translations
end

@testset "wrap_centers" begin
    wout = read_wout(joinpath(FIXTURE_PATH, "valence/band/silicon.wout"))

    centers = Wannier.wrap_centers(wout.centers, wout.lattice)

    ref_centers = [
        -3.35816 -4.72832 -4.71512 -3.40048
        3.35704 3.39463 4.71376 4.74103
        4.71664 3.36561 4.71388 3.40255
    ]
    ref_centers = [Vec3(ref_centers[:, i]) for i = 1:size(ref_centers,2)]

    @test all(isapprox.(centers, ref_centers; atol=2e-5))
end
