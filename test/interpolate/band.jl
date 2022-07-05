@testset "get_kpath_points" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    band = read_w90_bands(joinpath(FIXTURE_PATH, "valence/band/silicon"))

    # num points of 1st segment
    n_points = 100

    kpt_frac, x, symm_idx, symm_label = Wannier.get_kpath_points(
        win.kpoint_path,
        n_points,
        get_recip_lattice(win.unit_cell),
    )

    @test isapprox(kpt_frac, band.kpoints; atol = 1e-5)
    @test isapprox(x, band.x; atol = 1e-5)
    @test symm_idx == band.symm_idx
    @test symm_label == band.symm_label
end

@testset "interpolate" begin
    model = read_seedname(joinpath(FIXTURE_PATH, "valence/band/silicon"))
    band = read_w90_bands(joinpath(FIXTURE_PATH, "valence/band/silicon"))

    E = Wannier.interpolate(model, band.kpoints)

    @test isapprox(E, band.E; atol = 1e-8)
end
