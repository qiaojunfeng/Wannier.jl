# Shared setup for the SymmetryConstraint tests: Si2_hse has a 16-band window
# that cuts degenerate multiplets at several IBZ kpoints (nonsymmorphic group
# with time-reversal), which exercises the band-masking path of the projector.
@testitem "symmetry_constraint tables" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.symmetry_constraint(kstencil, isym, centers)

    # stars partition the FBZ mesh
    @test sort(vcat(sc.stars...)) == 1:sc.nk_fbz
    @test all(iki -> sc.fbz2ibz[sc.ibz2fbz[iki]][1] == iki, 1:sc.nk_ibz)

    # gauge expansion reproduces unfold_gauges
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.orbital_reps)
    f2i = get_kpoint_mappings(kstencil.kpoints, isym.kpoints_ibz, isym.symops)
    Ai = read_amn(dataset"Si2_hse/outputs/test/symmetrized.iamn").A
    Af_ref = Wannier.unfold_gauges(
        Ai, isym.kpoints_ibz, f2i, isym.symops, isym.orbital_reps, Rs
    )
    @test Wannier.expand_gauges(Ai, sc) ≈ Af_ref

    # cached overlap unfolding reproduces unfold_overlaps
    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    bvecs = get_bvectors(kstencil; fractional = true)
    Mf_ref, _, _ = Wannier.unfold_overlaps(
        mmn_i.M, mmn_i.kpb_k, mmn_i.kpb_G, isym.kpoints_ibz, bvecs,
        kstencil.kpoints, f2i, isym.spinors, isym.symops, isym.littlegroup_reps,
    )
    @test Wannier.unfold_overlaps_cached(mmn_i.M, sc) ≈ Mf_ref

    # pullback is the adjoint of the expansion: Re⟨G, expand(U)⟩ = Re⟨pullback(G), U⟩
    U = randn(ComplexF64, sc.nbands, sc.nwann, sc.nk_ibz)
    G = randn(ComplexF64, sc.nbands, sc.nwann, sc.nk_fbz)
    Gi = similar(U)
    Wannier.pullback_gauges!(Gi, G, sc)
    ip(a, b) = real(sum(conj.(a) .* b))
    @test isapprox(ip(G, Wannier.expand_gauges(U, sc)), ip(Gi, U); rtol = 1.0e-12)
end

@testitem "covariance projector" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.symmetry_constraint(kstencil, isym, centers)

    X = randn(ComplexF64, sc.nbands, sc.nwann, sc.nk_ibz)
    PX = Wannier.project_covariant(X, sc)
    PPX = Wannier.project_covariant(PX, sc)
    # Idempotence and single-average fixedness are limited by the numerical
    # quality of the d matrices in the isym file (~1e-5 for Si2_hse).
    @test norm(PPX - PX) < 1.0e-4
    @test norm(Wannier.covariant_average!(copy(PX), sc) - PX) < 1.0e-4

    # self-adjointness w.r.t. the real inner product (exact property)
    Y = randn(ComplexF64, size(X))
    PY = Wannier.project_covariant(Y, sc)
    ip(a, b) = real(sum(conj.(a) .* b))
    @test abs(ip(Y, PX) - ip(PY, X)) < 1.0e-10

    # covariant gauges carry no weight on symmetry-broken bands
    @test any(iki -> !all(sc.band_ok[iki]), 1:sc.nk_ibz)  # Si2_hse has some
    for iki in 1:sc.nk_ibz
        @test norm(PX[.!sc.band_ok[iki], :, iki]) < 1.0e-8
    end
end

@testitem "globalize_stencil" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    gs = Wannier.globalize_stencil(kstencil)
    bglob = get_bvectors(gs; fractional = true)
    @test bglob ≈ get_bvectors(kstencil; fractional = true)
    # every kpoint now uses the global b ordering
    for ik in 1:length(gs.kpoints), ib in 1:length(bglob)
        b = gs.kpoints[gs.kpb_k[ib, ik]] + gs.kpb_G[ib, ik] - gs.kpoints[ik]
        @test isapprox(b, bglob[ib]; atol = 1.0e-8)
    end
    @test gs.bweights ≈ kstencil.bweights
end
