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

    ref = readdlm(dataset"Si2_hse/outputs/test/kpt_map.txt", '\t', Int; comments=true)

    @test f2i == eachrow(ref)
end

@testitem "unfold_eigvals" begin
    using WannierIO
    using Wannier.Datasets

    win = read_win(dataset"Si2_hse/Si2.win")
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    kpoints_fbz = win.kpoints
    f2i = get_kpoint_mappings(kpoints_fbz, isym.kpoints_ibz, isym.symops)

    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    Ef = Wannier.unfold_eigvals(Ei, f2i)
    ref = read_eig(dataset"Si2_hse/Si2.eig")

    @test isapprox(Ef, ref)
end

@testitem "find_wf_symmetry_translations" begin
    using WannierIO, Wannier.Datasets
    using DelimitedFiles

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    isym = read_isym(dataset"Si2_hse/Si2.isym")

    centers = [p.center for p in nnkp.projections]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.repmat_wann)

    ref, header = readdlm(
        dataset"Si2_hse/outputs/test/R_translations.txt",
        ' ',
        Int;
        comments=true,
        header=true,
    )
    nsym, nwann = parse.(Int, header[1, 1:2])
    # Reshape ref into ref[isym][iwf][1:3]
    ref = [eachrow(ref[((is - 1) * nwann + 1):(is * nwann), :]) for is in 1:nsym]

    @test isapprox(Rs, ref)
end

@testitem "symmetrize_gauges" begin
    using WannierIO, Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    isym.repmat_band .= Wannier.rescale.(isym.repmat_band)

    centers = [p.center for p in nnkp.projections]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.repmat_wann)

    Ai = read_amn(dataset"Si2_hse/Si2.iamn")
    Asymm = Wannier.symmetrize_gauges(
        Ai, isym.kpoints_ibz, isym.symops, isym.repmat_band, isym.repmat_wann, Rs
    )

    ref = read_amn(dataset"Si2_hse/outputs/test/symmetrized.iamn")

    @test isapprox(Asymm, ref)
end

@testitem "unfold_gauges" begin
    using WannierIO, Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    isym.repmat_band .= Wannier.rescale.(isym.repmat_band)
    f2i = get_kpoint_mappings(kstencil.kpoints, isym.kpoints_ibz, isym.symops)

    centers = [p.center for p in nnkp.projections]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.repmat_wann)

    Asymm = read_amn(dataset"Si2_hse/outputs/test/symmetrized.iamn")
    Af = Wannier.unfold_gauges(
        Asymm, isym.kpoints_ibz, f2i, isym.symops, isym.repmat_wann, Rs
    )

    ref = read_amn(dataset"Si2_hse/Si2.amn")

    @test isapprox(Af, ref; atol=1e-10)
end
