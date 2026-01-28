@testitem "rescale" begin
    using WannierIO: RepMatBand

    d = [1.0 0.0; 0.0 1-1e-4]
    rep = RepMatBand{2}(1, 1, d)
    r = Wannier.rescale(rep)
    ref = [1.0 0.0; 0.0 1.0]

    @test isapprox(r.d, ref)
end

@testitem "find_wf_symmetry_translations" begin
    using WannierIO, Wannier.Datasets
    # sym = read_isym(dataset"Si2_hse/Si2.isym")
    sym = read_isym("/home/jqiao/git/WannierDatasets/datasets/Si2_hse/Si2.isym")

    # TODO
    @test isapprox(r, wout.centers; atol=1e-5)
end
