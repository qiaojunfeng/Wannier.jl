@testitem "rescale" begin
    using WannierIO: RepMatBand

    d = [1.0 0.0; 0.0 1 - 1.0e-4]
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
    kpoints_fbz = win["kpoints"]
    f2i = get_kpoint_mappings(kpoints_fbz, isym.kpoints_ibz, isym.symops)

    ref = readdlm(dataset"Si2_hse/outputs/test/kpt_map.txt", '\t', Int; comments = true)

    @test f2i == eachrow(ref)
end

@testitem "unfold_eigvals" begin
    using WannierIO
    using Wannier.Datasets

    win = read_win(dataset"Si2_hse/Si2.win")
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    kpoints_fbz = win["kpoints"]
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

    centers = [p.center for p in nnkp["projections"]]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.repmat_wann)

    ref, header = readdlm(
        dataset"Si2_hse/outputs/test/R_translations.txt",
        ' ',
        Int;
        comments = true,
        header = true,
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
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.repmat_band)

    centers = [p.center for p in nnkp["projections"]]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.repmat_wann)

    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    Asymm = Wannier.symmetrize_gauges(
        Ai, isym.kpoints_ibz, isym.symops, isym.repmat_band, isym.repmat_wann, Rs
    )

    ref = read_amn(dataset"Si2_hse/outputs/test/symmetrized.iamn").A

    @test isapprox(Asymm, ref)
end

@testitem "unfold_gauges" begin
    using WannierIO, Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.repmat_band)
    f2i = get_kpoint_mappings(kstencil.kpoints, isym.kpoints_ibz, isym.symops)

    centers = [p.center for p in nnkp["projections"]]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.repmat_wann)

    Asymm = read_amn(dataset"Si2_hse/outputs/test/symmetrized.iamn").A
    Af = Wannier.unfold_gauges(
        Asymm, isym.kpoints_ibz, f2i, isym.symops, isym.repmat_wann, Rs
    )

    ref = read_amn(dataset"Si2_hse/Si2.amn").A

    @test isapprox(Af, ref; atol = 1.0e-10)
end

@testitem "get_equivalence_mappings" begin
    using WannierIO, Wannier.Datasets
    using DelimitedFiles

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    symops = read_isym(dataset"Si2_hse/Si2.isym").symops
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    bvecs = get_bvectors(kstencil; fractional = true)

    equiv = Wannier.get_equivalence_mappings(bvecs, symops)
    ref = readdlm(dataset"Si2_hse/outputs/test/b_equivalence.txt", ' ', Int)
    @test equiv == ref
end

@testitem "merge_symops" begin
    using WannierIO, Wannier.Datasets

    isym = read_isym(dataset"Si2_hse/Si2.isym")
    isym_kbi, isym_kf, isym_kbf = 6, 11, 2

    isym_h, factor, T = Wannier.merge_symops(
        isym.spinors, isym.symops, [isym_kbi, isym_kf, isym_kbf], [true, true, false]
    )

    @test isym_h == 24
    @test factor == 1
    @test T == [0, 0, 0]
end

@testitem "unfold_overlaps" begin
    using LinearAlgebra
    using WannierIO, Wannier, Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.repmat_band)
    f2i = get_kpoint_mappings(kstencil.kpoints, isym.kpoints_ibz, isym.symops)

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Mi = mmn_i.M
    kpb_k_i = mmn_i.kpb_k
    kpb_G_i = mmn_i.kpb_G

    bvectors = get_bvectors(kstencil; fractional = true)
    Mf, kpb_k_f, kpb_G_f = Wannier.unfold_overlaps(
        Mi,
        kpb_k_i,
        kpb_G_i,
        isym.kpoints_ibz,
        bvectors,
        kstencil.kpoints,
        f2i,
        isym.spinors,
        isym.symops,
        isym.repmat_band,
    )

    mmn_ref = read_mmn(dataset"Si2_hse/Si2.mmn")
    Mref = mmn_ref.M
    kpb_k_ref = mmn_ref.kpb_k
    kpb_G_ref = mmn_ref.kpb_G

    # There is approx 1e-8 differences compared to the reference mmn, because
    # the reference mmn was generated with python code, where the b vectors
    # has slight numerical errors (should be e.g. 0.125 (as in julia) but in
    # python it was 0.1249999962372273), therefore the phase factor was slightly
    # wrong.
    d = map(Mf - Mref) do mk
        map(mk) do mkb
            maximum(norm, mkb)
        end
    end
    dmax = maximum(maximum.(d))
    @test all(isapprox(0; atol = 1.0e-7), dmax)
    @test kpb_k_f == kpb_k_ref
    @test kpb_G_f == kpb_G_ref
end
