
@testset "Rvectors WS" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    ref_Rvecs, H, positions = Wannier.read_w90_tb(
        joinpath(FIXTURE_PATH, "valence/band/ws/silicon")
    )

    Rvecs = Wannier.get_Rvectors_ws(ref_Rvecs.lattice, win.mp_grid)

    @test Rvecs.R == ref_Rvecs.R
    @test Rvecs.N == ref_Rvecs.N
end

@testset "Rvectors MDRS" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    ref_Rvecs, H, positions = Wannier.read_w90_tb(
        joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon")
    )
    wout = read_wout(joinpath(FIXTURE_PATH, "valence/band/silicon.wout"))

    lattice = ref_Rvecs.lattice
    # to fractional coordinates
    centers = inv(lattice) * wout.centers
    Rvecs = Wannier.get_Rvectors_mdrs(lattice, win.mp_grid, centers)

    @test Rvecs.R == ref_Rvecs.R
    @test Rvecs.N == ref_Rvecs.N
    @test Rvecs.T == ref_Rvecs.T
    @test Rvecs.Nᵀ == ref_Rvecs.Nᵀ
end
