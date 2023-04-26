
@testset "write chk from Model" begin
    chk = WannierIO.read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))
    model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))
    model.U .= WannierIO.get_U(chk)

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, model)

    chk2 = WannierIO.read_chk(tmpfile)

    # header contains date, always different
    # @test chk.header == chk2.header
    @test chk.exclude_bands == chk2.exclude_bands
    @test chk.lattice ≈ chk2.lattice
    @test chk.recip_lattice ≈ chk2.recip_lattice
    @test chk.kgrid == chk2.kgrid
    @test chk.kpoints ≈ chk2.kpoints
    @test chk.checkpoint == chk2.checkpoint
    @test chk.have_disentangled == chk2.have_disentangled
    @test chk.ΩI ≈ chk2.ΩI
    # the dis_bands are all true when writing from a model
    @test chk2.dis_bands == trues(chk2.n_bands, chk2.n_kpts)
    @test model.U ≈ WannierIO.get_U(chk2)

    @test chk.M ≈ chk2.M
    @test chk.r ≈ chk2.r
    @test chk.ω ≈ chk2.ω
end

@testset "Model(chk)" begin
    chk = WannierIO.read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))
    model = Wannier.Model(chk)

    @test chk.lattice ≈ model.lattice
    @test chk.recip_lattice ≈ model.recip_lattice
    @test chk.kgrid == model.kgrid
    @test chk.kpoints ≈ model.kpoints
    @test chk.M ≈ model.M
end
