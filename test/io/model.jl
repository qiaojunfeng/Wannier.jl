@testitem "read_w90" begin
    using Wannier.Datasets
    model = load_dataset("Si2")

    @test n_bands(model) == 16
    @test n_wannier(model) == 8
    @test n_kpoints(model) == 9^3
end

@testitem "write_w90" begin
    using Wannier.Datasets
    model = load_dataset("Si2")

    outdir = mktempdir(; cleanup=true)
    outprefix = joinpath(outdir, "silicon")
    write_w90(outprefix, model)

    gauges = read_amn_ortho("$outprefix.amn")
    @test gauges ≈ model.gauges

    eigenvalues = read_eig("$outprefix.eig")
    @test eigenvalues ≈ model.eigenvalues

    overlaps, kpb_k, kpb_G = read_mmn("$outprefix.mmn")
    @test overlaps ≈ model.overlaps
    @test kpb_k == model.kpb_k
    @test kpb_G == model.kpb_G
end

@testitem "read_w90_with_chk" begin
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/reference/Si2.chk")

    @test model isa Wannier.Model
    @test n_wannier(model) == 8
end

@testitem "write chk from Model" begin
    using LinearAlgebra
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/reference/Si2.chk.fmt")
    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, model)

    chk = read_chk(tmpfile)

    # header contains date, always different
    # @test chk.header == chk2.header
    @test isempty(chk.exclude_bands)
    @test model.lattice ≈ chk.lattice
    @test model.recip_lattice ≈ chk.recip_lattice
    @test model.kgrid_size == chk.kgrid
    @test model.kpoints ≈ chk.kpoints
    @test true == chk.have_disentangled
    @test model.entangled_bands == chk.dis_bands
    # the Hamiltonian rotated by Udis must be diagonal, according to W90 convention
    H = transform_gauge(model.eigenvalues, Wannier.get_Udis(chk))
    # this is too strict, even
    #   norm(H[:, :, ik] - Hdiag[:, :, ik]) ≈ 1e-14
    # is still false
    # @test all(isdiag(H[:, :, ik]) for ik in axes(H, 3))
    Hdiag = map(H) do h
        Hermitian(Diagonal(h))
    end
    @test H ≈ Hdiag
    # the unitary matrix should be the same
    @test model.gauges ≈ Wannier.get_U(chk)

    M = transform_gauge(model.overlaps, model.kpb_k, model.gauges)
    @test M ≈ chk.M
    Ω = omega(model, model.gauges)
    @test Ω.ΩI ≈ chk.ΩI
    @test Ω.r ≈ chk.r
    @test Ω.ω ≈ chk.ω
end

@testitem "Model from chk" begin
    using Wannier.Datasets
    chk = read_chk(dataset"Si2/reference/Si2.chk")
    model = Wannier.Model(chk)

    @test chk.lattice ≈ model.lattice
    @test chk.recip_lattice ≈ reciprocal_lattice(model)
    @test chk.kgrid == model.kgrid_size
    @test chk.kpoints ≈ model.kpoints
    @test chk.M ≈ model.overlaps
end
