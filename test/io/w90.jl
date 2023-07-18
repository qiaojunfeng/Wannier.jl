using YAML
using Brillouin

@testset "read nnkp" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/nnkp.yaml")

    bvectors = read_nnkp(joinpath(FIXTURE_PATH, "silicon/silicon.nnkp"))

    # Convert type so YAML can write it.
    kpb_G = bvectors.kpb_G
    dict = Dict(
        "recip_lattice" => [[bvectors.recip_lattice[i, j] for i in 1:3] for j in 1:3],
        "kpoints" => bvectors.kpoints,
        "bvectors" => bvectors.bvectors,
        "kpb_k" => bvectors.kpb_k,
        "kpb_G" => bvectors.kpb_G,
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
    @test bvectors.kpb_G ≈ bvectors2.kpb_G
end

@testset "read/write w90 band" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    recip_lattice = Wannier.get_recip_lattice(win.unit_cell_cart)
    kpi, eigenvalues = read_w90_band(
        joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"), recip_lattice
    )

    outdir = mktempdir(; cleanup=true)
    outseedname = joinpath(outdir, "silicon")

    write_w90_band(outseedname, kpi, eigenvalues)

    kpi2, eigenvalues2 = read_w90_band(outseedname, recip_lattice)

    @test kpi.kpaths ≈ kpi2.kpaths
    @test kpi.labels == kpi2.labels
    @test kpi.basis ≈ kpi2.basis
    @test Symbol(kpi.setting) == Symbol(kpi2.setting)
    @test eigenvalues ≈ eigenvalues2
end

@testset "read tb" begin
    tb_ws = Wannier.read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/ws/silicon"))
    Rvecs = tb_ws.Rvectors
    # just some simple tests
    R1 = [-3, 1, 1]
    @test Rvecs.R[1] == R1 == tb_ws.H[1].R_cryst
    H111 = 0.51893360E-02 + im * -0.29716277E-02
    @test tb_ws.H[1][1, 1] ≈ H111
    P111end = 0.24832468E-03 + im * -0.21054981E-03
    @test tb_ws.r_x[end][1, 1] ≈ P111end

    tb_mdrs = Wannier.read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))
    Rvecs = tb_mdrs.Rvectors
    @test Rvecs.R[1] == R1
    @test Rvecs.T[1][1, 1] == [Vec3(0, 0, 0), Vec3(4, -4, 0), Vec3(4, 0, -4), Vec3(4, 0, 0)]
    @test Rvecs.Nᵀ[1][1, 1] == 4
    H111 = 0.0012973340069440233 - 0.0007429069229594317im
    @test tb_mdrs.H[1][1, 1] ≈ H111
    P111end = 0.0 + 0.0im
    @test tb_mdrs.r_x[end][1, 1] ≈ P111end
end

@testset "read tb kpoints/kpath" begin
    win = Wannier.read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    interp_model = Wannier.read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon");)

    @test interp_model.H isa Wannier.TBHamiltonian
end
