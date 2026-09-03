@testitem "rescale" begin
    using WannierIO: LittleGroupRep

    d = [1.0 0.0; 0.0 1 - 1.0e-4]
    rep = LittleGroupRep{2}(1, 1, d)
    r = Wannier.rescale(rep)
    ref = [1.0 0.0; 0.0 1.0]

    @test isapprox(r.d, ref)
end

@testitem "map_fbz_to_ibz" begin
    using WannierIO
    using Wannier.Datasets
    using DelimitedFiles

    win = read_win(dataset"Si2_hse/Si2.win")
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    kpoints_fbz = win["kpoints"]
    f2i = map_fbz_to_ibz(kpoints_fbz, isym.kpoints_ibz, isym.symops)

    ref = readdlm(dataset"Si2_hse/outputs/test/kpt_map.txt", '\t', Int; comments = true)

    @test f2i == eachrow(ref)
end

@testitem "unfold_eigvals" begin
    using WannierIO
    using Wannier.Datasets

    win = read_win(dataset"Si2_hse/Si2.win")
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    kpoints_fbz = win["kpoints"]
    f2i = map_fbz_to_ibz(kpoints_fbz, isym.kpoints_ibz, isym.symops)

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
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.orbital_reps)

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

    # The reference file stores the translations of the *exact inverse*
    # operations indexed by the forward operation (the historical
    # convention): ref[is] = R(g_is^{-1}). In the standard convention
    # `Rs[j] = R(g_j)` of the operation itself, related by
    # Rs[invs(is)] = ref[is] + L(is), where g_is^{-1} = t_{-L} ∘ g_{invs(is)}.
    for is in 1:nsym
        isinv = isym.symops[is].isym_inv
        L = Wannier.inverse_translation_mismatch(isym.symops, is)
        @test all(Rs[isinv][iw] == ref[is][iw] + L for iw in 1:nwann)
    end
end

@testitem "symmetrize_gauges" begin
    using WannierIO, Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)

    centers = [p.center for p in nnkp["projections"]]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.orbital_reps)

    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    Asymm = Wannier.symmetrize_gauges(
        Ai, isym.kpoints_ibz, isym.symops, isym.littlegroup_reps, isym.orbital_reps, Rs
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
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    f2i = map_fbz_to_ibz(kstencil.kpoints, isym.kpoints_ibz, isym.symops)

    centers = [p.center for p in nnkp["projections"]]
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.orbital_reps)

    Asymm = read_amn(dataset"Si2_hse/outputs/test/symmetrized.iamn").A
    Af = Wannier.unfold_gauges(
        Asymm, isym.kpoints_ibz, f2i, isym.symops, isym.orbital_reps, Rs
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

@testitem "compose_symops" begin
    using WannierIO, Wannier.Datasets

    isym = read_isym(dataset"Si2_hse/Si2.isym")
    isym_kbi, isym_kf, isym_kbf = 6, 11, 2

    isym_h, factor, T = Wannier.compose_symops(
        isym.spinors, isym.symops, [isym_kbi, isym_kf, isym_kbf], [true, true, false]
    )

    @test isym_h == 24
    @test factor == 1
    @test T == [0, 0, 0]
end

@testitem "reconstruct_overlaps" begin
    using LinearAlgebra
    using WannierIO, Wannier, Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    f2i = map_fbz_to_ibz(kstencil.kpoints, isym.kpoints_ibz, isym.symops)

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Mi = mmn_i.M
    kpb_k_i = mmn_i.kpb_k
    kpb_G_i = mmn_i.kpb_G

    bvectors = get_bvectors(kstencil; fractional = true)
    Mf, kpb_k_f, kpb_G_f = Wannier.reconstruct_overlaps(
        Mi,
        kpb_k_i,
        kpb_G_i,
        isym.kpoints_ibz,
        bvectors,
        kstencil.kpoints,
        f2i,
        isym.spinors,
        isym.symops,
        isym.littlegroup_reps,
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
