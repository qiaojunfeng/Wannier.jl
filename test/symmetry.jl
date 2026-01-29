@testitem "rescale" begin
    using WannierIO: RepMatBand

    d = [1.0 0.0; 0.0 1-1e-4]
    rep = RepMatBand{2}(1, 1, d)
    r = Wannier.rescale(rep)
    ref = [1.0 0.0; 0.0 1.0]

    @test isapprox(r.d, ref)
end

@testitem "get_kpoint_mappings" begin
    using WannierIO
    using Wannier.Datasets
    using DelimitedFiles

    win = read_win(dataset"Si2_hse/Si2.win")
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    kpoints_fbz = win.kpoints
    f2i = get_kpoint_mappings(kpoints_fbz, isym.kpoints_ibz, isym.symops)

    ref = readdlm(dataset"Si2_hse/outputs/kpt_map.txt", '\t', Int; comments=true)

    @test f2i == eachrow(ref)
end

@testitem "find_wf_symmetry_translations" begin
    using WannierIO, Wannier.Datasets
    # sym = read_isym(dataset"Si2_hse/Si2.isym")
    sym = read_isym("/home/jqiao/git/WannierDatasets/datasets/Si2_hse/Si2.isym")

    # TODO
    @test isapprox(r, wout.centers; atol=1e-5)
end
