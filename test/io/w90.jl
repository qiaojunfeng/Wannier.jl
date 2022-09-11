using YAML
using Brillouin

@testset "read nnkp" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/nnkp.yaml")

    bvectors = read_nnkp(joinpath(FIXTURE_PATH, "silicon/silicon.nnkp"))

    # Convert type so YAML can write it.
    kpb_b = bvectors.kpb_b
    dict = Dict(
        "recip_lattice" => mat2vec(bvectors.recip_lattice),
        "kpoints" => mat2vec(bvectors.kpoints),
        "bvectors" => mat2vec(bvectors.bvectors),
        "kpb_k" => mat2vec(bvectors.kpb_k),
        "kpb_b" => [mat2vec(kpb_b[:, :, ik]) for ik in axes(kpb_b, 3)],
    )

    # YAML.write_file(String(@__DIR__) * "/test_data/nnkp.yaml", dict)

    for (key, value) in dict
        @test value ≈ test_data[key]
    end
end

@testset "read/write nnkp" begin
    bvectors = read_nnkp(joinpath(FIXTURE_PATH, "silicon/silicon.nnkp"))
    tmpfile = tempname(; cleanup=true)
    n_wann = 8
    write_nnkp(tmpfile, bvectors, n_wann)

    bvectors2 = read_nnkp(tmpfile)
    @test bvectors.recip_lattice ≈ bvectors2.recip_lattice
    @test bvectors.kpoints ≈ bvectors2.kpoints
    @test bvectors.bvectors ≈ bvectors2.bvectors
    @test bvectors.kpb_k ≈ bvectors2.kpb_k
    @test bvectors.kpb_b ≈ bvectors2.kpb_b
end

@testset "read/write w90 band" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    recip_lattice = Wannier.get_recip_lattice(win.unit_cell)
    kpi, E = read_w90_band(
        joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"), recip_lattice
    )

    outdir = mktempdir(; cleanup=true)
    outseedname = joinpath(outdir, "silicon")

    write_w90_band(outseedname, kpi, E)

    kpi2, E2 = read_w90_band(outseedname, recip_lattice)

    @test kpi.kpaths ≈ kpi2.kpaths
    @test kpi.labels == kpi2.labels
    @test kpi.basis ≈ kpi2.basis
    @test Symbol(kpi.setting) == Symbol(kpi2.setting)
    @test E ≈ E2
end

@testset "read tb" begin
    Rvecs, H, positions = Wannier.read_w90_tb(
        joinpath(FIXTURE_PATH, "valence/band/ws/silicon")
    )
    # just some simple tests
    R1 = [-3, 1, 1]
    @test Rvecs.R[:, 1] == R1
    H111 = 0.51893360E-02 + im * -0.29716277E-02
    @test H[1, 1, 1] ≈ H111
    P111end = 0.24832468E-03 + im * -0.21054981E-03
    @test positions[1, 1, 1, end] ≈ P111end

    Rvecs, H, positions = Wannier.read_w90_tb(
        joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon")
    )
    @test Rvecs.R[:, 1] == R1
    @test Rvecs.T[1, 1, 1] == [0 4 4 4; 0 -4 0 0; 0 0 -4 0]
    @test Rvecs.Nᵀ[1, 1, 1] == 4
    @test H[1, 1, 1] ≈ H111
    @test positions[1, 1, 1, end] ≈ P111end
end
