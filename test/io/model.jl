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
    @test kpb_k == model.bvectors.kpb_k
    @test kpb_G == model.bvectors.kpb_G
end

@testitem "read_w90_with_chk" begin
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/reference/Si2.chk")

    @test model isa Wannier.Model
    @test n_wannier(model) == 8
end
