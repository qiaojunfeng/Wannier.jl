@testitem "find_wf_symmetry_translations" begin
    using WannierIO, Wannier.Datasets
    # sym = read_isym(dataset"Si2_hse/Si2.isym")
    sym = read_isym("/home/jqiao/git/WannierDatasets/datasets/Si2_hse/Si2.isym")

    # TODO
    @test isapprox(r, wout.centers; atol=1e-5)
end
