# Shared setup for the SymmetryConstraint tests: Si2_hse has a 16-band window
# that cuts degenerate multiplets at several IBZ kpoints (nonsymmorphic group
# with time-reversal), which exercises the band-masking path of the projector.
@testitem "SymmetryConstraint tables" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(kstencil, isym, centers)

    # stars partition the FBZ mesh
    @test sort(vcat(sc.stars...)) == 1:sc.nk_fbz
    @test all(iki -> sc.fbz2ibz[sc.ibz2fbz[iki]][1] == iki, 1:sc.nk_ibz)

    # gauge expansion reproduces unfold_gauges
    Rs = Wannier.find_wf_symmetry_translations(centers, isym.symops, isym.orbital_reps)
    f2i = map_fbz_to_ibz(kstencil.kpoints, isym.kpoints_ibz, isym.symops)
    Ai = read_amn(dataset"Si2_hse/outputs/test/symmetrized.iamn").A
    Af_ref = Wannier.unfold_gauges(
        Ai, isym.kpoints_ibz, f2i, isym.symops, isym.orbital_reps, Rs
    )
    @test Wannier.reconstruct_gauges(Ai, sc) ≈ Af_ref

    # cached overlap reconstruction reproduces reconstruct_overlaps
    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    bvecs = get_bvectors(kstencil; fractional = true)
    Mf_ref, _, _ = Wannier.reconstruct_overlaps(
        mmn_i.M, mmn_i.kpb_k, mmn_i.kpb_G, isym.kpoints_ibz, bvecs,
        kstencil.kpoints, f2i, isym.spinors, isym.symops, isym.littlegroup_reps,
    )
    @test Wannier.reconstruct_overlaps(mmn_i.M, sc) ≈ Mf_ref

    # pullback is the adjoint of the expansion: Re⟨G, expand(U)⟩ = Re⟨pullback(G), U⟩
    U = randn(ComplexF64, sc.nbands, sc.nwann, sc.nk_ibz)
    G = randn(ComplexF64, sc.nbands, sc.nwann, sc.nk_fbz)
    Gi = similar(U)
    Wannier.pullback_gauges!(Gi, G, sc)
    ip(a, b) = real(sum(conj.(a) .* b))
    @test isapprox(ip(G, Wannier.reconstruct_gauges(U, sc)), ip(Gi, U); rtol = 1.0e-12)
end

@testitem "covariance projector" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    clean_littlegroup_reps!(isym.littlegroup_reps, Ei)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(kstencil, isym, centers)

    X = randn(ComplexF64, sc.nbands, sc.nwann, sc.nk_ibz)
    PX = Wannier.project_covariant(X, sc)
    PPX = Wannier.project_covariant(PX, sc)
    # Idempotence and single-average fixedness are limited by the numerical
    # quality of the d matrices (raw Si2_hse carries ~1e-5 noise; cleaning
    # lowers the floor to ~1e-7, limited by the masked broken multiplets).
    @test norm(PPX - PX) < 1.0e-5
    @test norm(Wannier.group_average!(copy(PX), sc) - PX) < 1.0e-5

    # self-adjointness w.r.t. the real inner product (exact property)
    Y = randn(ComplexF64, size(X))
    PY = Wannier.project_covariant(Y, sc)
    ip(a, b) = real(sum(conj.(a) .* b))
    @test abs(ip(Y, PX) - ip(PY, X)) < 1.0e-10

    # covariant gauges carry no weight on symmetry-broken bands
    @test any(iki -> !all(sc.covariant_bands[iki]), 1:sc.nk_ibz)  # Si2_hse has some
    for iki in 1:sc.nk_ibz
        @test norm(PX[.!sc.covariant_bands[iki], :, iki]) < 1.0e-8
    end
end

@testitem "globalize_bvector_ordering" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    gs = Wannier.globalize_bvector_ordering(kstencil)
    bglob = get_bvectors(gs; fractional = true)
    @test bglob ≈ get_bvectors(kstencil; fractional = true)
    # every kpoint now uses the global b ordering
    for ik in 1:length(gs.kpoints), ib in 1:length(bglob)
        b = gs.kpoints[gs.kpb_k[ib, ik]] + gs.kpb_G[ib, ik] - gs.kpoints[ik]
        @test isapprox(b, bglob[ib]; atol = 1.0e-8)
    end
    @test gs.bweights ≈ kstencil.bweights
end

@testitem "opposite b vectors on a singleton mesh axis" begin
    # Along a singleton mesh axis, +G and -G are equivalent k-point shifts but
    # remain distinct displacement vectors in the finite-difference shell.
    bvecs = [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1 / 7, 0.0, 0.0],
        [-1 / 7, 0.0, 0.0],
    ]
    minus_b = Wannier._opposite_bvector_indices(bvecs)
    @test minus_b == [2, 1, 4, 3]
    @test minus_b[minus_b] == eachindex(bvecs)
end

@testitem "spinor transport sign is carried once" begin
    d = ComplexF64[1 + 2im 3 + 4im; 5 + 6im 7 + 8im]
    factor = -1
    θ1, θ2 = 1 / 8, -1 / 4
    phase_ref = factor * exp(-im * 2π * (θ1 + θ2))

    phase, dmat = Wannier._overlap_transport_seed(factor, θ1, θ2, d, false)
    @test phase ≈ phase_ref
    @test dmat == d
    @test phase .* dmat ≈ phase_ref .* d

    phase_trev, dmat_trev = Wannier._overlap_transport_seed(factor, θ1, θ2, d, true)
    @test phase_trev ≈ phase_ref
    @test dmat_trev == conj.(d)
    @test phase_trev .* dmat_trev ≈ phase_ref .* conj.(d)
end

@testitem "fullmesh symmetric fg" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(ks0, isym, centers)
    ks = Wannier.globalize_bvector_ordering(ks0)

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Mf = Wannier.reconstruct_overlaps(mmn_i.M, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    atom_positions = [p.second for p in win["atoms_frac"]]
    atom_labels = map(x -> string(x.first), win["atoms_frac"])
    model = Wannier.Model(
        win["unit_cell_cart"], atom_positions, atom_labels,
        ks, Mf, Wannier.reconstruct_gauges(Wannier.project_covariant(Ai, sc), sc), Ef, frozen,
    )

    sm = Wannier.SymmetricModel(model, sc, mmn_i.M)
    prob = Wannier.Problem(
        Wannier.Variance(), sm, Wannier.SymmetricXYLayout(:fullmesh)
    )
    fg! = Wannier._optimizer_callback(prob)
    x = Wannier.initial_parameters(prob.layout, sm)
    G = zero(x)
    Ω = fg!(1.0, G, x)

    # value consistency with the plain full-mesh spread of the reconstructed gauge
    Ufull, _ = Wannier.finalize_result(prob.layout, x, sm)
    @test Ω ≈ Wannier.spread(model.kstencil, model.overlaps, Ufull).Ω

    # Directional finite-difference gradient check. The compact XY layout does
    # not store the entries fixed by the frozen subspace. The rtol is limited
    # by FD truncation of the Si2_hse point (near-branch-point
    # Im-log diagonals); the same check on clean data (Ge4Ru4) passes at 1e-9.
    f = x -> fg!(1.0, nothing, x)
    for _ in 1:2
        dx = randn(ComplexF64, size(x))
        dx ./= norm(dx)
        ε = 1.0e-4
        fd = (f(x .+ ε .* dx) - f(x .- ε .* dx)) / (2ε)
        an = real(sum(conj.(G) .* dx))
        @test isapprox(fd, an; rtol = 1.0e-4)
    end
end

@testitem "transport equals fullmesh" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    clean_littlegroup_reps!(isym.littlegroup_reps, Ei)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(ks0, isym, centers)
    ks = Wannier.globalize_bvector_ordering(ks0)

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Mf = Wannier.reconstruct_overlaps(mmn_i.M, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    atom_positions = [p.second for p in win["atoms_frac"]]
    atom_labels = map(x -> string(x.first), win["atoms_frac"])
    model = Wannier.Model(
        win["unit_cell_cart"], atom_positions, atom_labels,
        ks, Mf, Wannier.reconstruct_gauges(Wannier.project_covariant(Ai, sc), sc), Ef, frozen,
    )

    sm = Wannier.SymmetricModel(model, sc, mmn_i.M)
    prob1 = Wannier.Problem(
        Wannier.Variance(), sm, Wannier.SymmetricXYLayout(:fullmesh)
    )
    prob2 = Wannier.Problem(
        Wannier.Variance(), sm, Wannier.SymmetricXYLayout(:transport)
    )
    fg1! = Wannier._optimizer_callback(prob1)
    fg2! = Wannier._optimizer_callback(prob2)
    x = Wannier.initial_parameters(prob2.layout, sm)
    G1, G2 = zero(x), zero(x)
    Ω1 = fg1!(1.0, G1, x)
    Ω2 = fg2!(1.0, G2, x)

    # The transport theorem is exact for a covariant gauge; the agreement is
    # limited only by the projector's data-precision fixed point. With
    # `clean_littlegroup_reps!` the d matrices are exact (raw Si2_hse noise
    # ~1e-5 would floor the agreement at ~1e-7 / ~1e-8).
    @test isapprox(Ω1, Ω2; atol = 1.0e-10)
    @test norm(G1 - G2) / norm(G1) < 1.0e-11

    # transport identity for the Wannier-gauge overlaps themselves:
    # M̃(kf, bf) = phase · K_f[L† M̃_i R] versus the full-mesh product
    U_fbz, U_ibz = Wannier.finalize_result(prob2.layout, x, sm)
    ikf = sc.stars[2][end]   # a star member away from its IBZ representative
    iki = sc.fbz2ibz[ikf][1]
    for ibf in 1:sc.nbvecs
        ikpb = ks.kpb_k[ibf, ikf]
        Mt_direct = U_fbz[:, :, ikf]' * Mf[:, :, ibf, ikf] * U_fbz[:, :, ikpb]
        ibi = sc.ibi_of[ibf, ikf]
        ikb = sc.ikb[ibi, iki]
        Ukb = Wannier._kconj(U_ibz[:, :, ikb] * sc.Aib[ibi, iki], sc.trev_kb[ibi, iki])
        Mt_i = U_ibz[:, :, iki]' * mmn_i.M[:, :, ibi, iki] * Ukb
        Mt_trans = sc.phase[ibf, ikf] .*
            Wannier._kconj(Matrix(sc.Lmat[ikf])' * Mt_i * Matrix(sc.Rmat[ibf, ikf]), sc.trev_f[ikf])
        @test norm(Mt_trans - Mt_direct) < 1.0e-5
    end
end

@testitem "hermiticity-pair transport" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(ks0, isym, centers)

    # table consistency: minus_b is the −b involution, ikpb_fbz points at ki+bi,
    # and the dagger member star-maps to the pass-0 IBZ point kb
    bvecs = get_bvectors(ks0; fractional = true)
    @test all(ib -> Wannier.isequivalent(bvecs[sc.minus_b[ib]], -bvecs[ib]), 1:sc.nbvecs)
    @test sc.minus_b[sc.minus_b] == 1:sc.nbvecs
    for iki in 1:sc.nk_ibz, ibi in 1:sc.nbvecs
        @test Wannier.isequivalent(
            ks0.kpoints[sc.ikpb_fbz[ibi, iki]], isym.kpoints_ibz[iki] + bvecs[ibi]
        )
        @test sc.fbz2ibz[sc.ikpb_fbz[ibi, iki]][1] == sc.ikb[ibi, iki]
    end

    # partner-transport identity at a covariant gauge: with p = P(q) the
    # Hermiticity pair of q = (iki, ibi),
    #   M̃_q = conj(phase2) · 𝒦_2[ R2† · M̃_p† · L2 ]
    # (the relation `_fg_transport_core!` uses in its loop C). For the 2-cycle pairs
    # of P — the only ones the halving derives — the stored overlaps are
    # Hermiticity-consistent and the identity is machine-exact even on this
    # noisy dataset (max ~5e-12; tested at 1e-4 for margin). Along P-chains
    # the two data entries are related only through a little-group rotation
    # and the residual is the data's symmetry noise (up to ~7e-4 here, ~4e-7
    # on clean Ge4Ru4 data); those pairs are computed directly.
    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    U_i = Wannier.project_covariant(Wannier.lowdin_orthonormalize(Ai), sc)
    nw, nbv, nki = sc.nwann, sc.nbvecs, sc.nk_ibz
    Mt = zeros(ComplexF64, nw, nw, nbv, nki)
    for iki in 1:nki, ibi in 1:nbv
        kb = sc.ikb[ibi, iki]
        Ukb = Wannier._kconj(U_i[:, :, kb] * sc.Aib[ibi, iki], sc.trev_kb[ibi, iki])
        Mt[:, :, ibi, iki] = U_i[:, :, iki]' * mmn_i.M[:, :, ibi, iki] * Ukb
    end
    partner(iki, ibi) = (sc.ikb[ibi, iki], sc.ibi_of[sc.minus_b[ibi], sc.ikpb_fbz[ibi, iki]])
    n2cycle = 0
    for iki in 1:nki, ibi in 1:nbv
        kb, ibi2 = partner(iki, ibi)
        partner(kb, ibi2) == (iki, ibi) || continue  # 2-cycle pairs only
        global n2cycle += 1
        ikpb = sc.ikpb_fbz[ibi, iki]
        ib2 = sc.minus_b[ibi]
        pred = conj(sc.phase[ib2, ikpb]) .* Wannier._kconj(
            Matrix(sc.Rmat[ib2, ikpb])' * Mt[:, :, ibi2, kb]' * Matrix(sc.Lmat[ikpb]),
            sc.trev_f[ikpb],
        )
        @test norm(pred - Mt[:, :, ibi, iki]) < 1.0e-4
    end
    # the halving must actually kick in: Si2_hse has 144 2-cycle pairs of 232
    @test n2cycle > 0.5 * nbv * nki
end

@testitem "schur block parametrization" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(ks0, isym, centers)

    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    frozen_ibz = frozen[:, sc.ibz2fbz]

    sb = Wannier.schur_basis(sc, frozen_ibz)
    # parameter reduction and feasibility (schur_basis errors when infeasible);
    # nx counts REAL parameters, the XY reference is complex
    @test 0 < sb.nx < 2 * (sc.nwann^2 + sc.nbands * sc.nwann) * sc.nk_ibz

    # anti-unitary metadata: pairing partner links are symmetric, derived
    # blocks carry no parameters, block shapes match across a pair, and
    # quaternionic blocks have even multiplicities
    for (iki, blks) in enumerate(sb.blocks), (ic, b) in enumerate(blks)
        if b.akind == Wannier.ANTIUNITARY_PAIRING_SOURCE
            p = blks[b.partner]
            @test p.akind == Wannier.ANTIUNITARY_PAIRING_DERIVED && p.partner == ic
            @test (p.dim, p.mb, p.mo, p.mf) == (b.dim, b.mb, b.mo, b.mf)
        elseif b.akind != Wannier.ANTIUNITARY_PAIRING_DERIVED
            @test b.partner == 0
        end
        if b.akind == Wannier.ANTIUNITARY_WIGNER_QUATERNIONIC
            @test iseven(b.mf) && iseven(b.mo) && iseven(b.mb)
        end
    end
    @test sb.nx == sum(
        b.akind == Wannier.ANTIUNITARY_PAIRING_DERIVED ? 0 :
            (
                b.akind == Wannier.ANTIUNITARY_WIGNER_REAL ||
                b.akind == Wannier.ANTIUNITARY_WIGNER_QUATERNIONIC ? 1 : 2
            ) *
            (b.mo^2 + (b.mb - b.mf) * (b.mo - b.mf))
            for blks in sb.blocks for b in blks
    )

    # assembly of the initial parameters: covariant (to the isym data noise),
    # semi-unitary, and consistent with the transport-path value at the same gauge
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    U0 = Wannier.project_covariant(Ai, sc)
    x0 = Wannier.initial_parameters(U0, sb)
    @test eltype(x0) <: Real   # re/im pairs; real/quaternion components
    Ud = zeros(ComplexF64, sc.nbands, sc.nwann, sc.nk_ibz)
    Wannier.assemble_gauge!(Ud, x0, sb)
    @test Wannier.covariance_residual(Ud, sc) < 1.0e-3
    # semi-unitarity of the assembly is limited by the isym data noise through
    # the anti-unitary coset average (1.9e-12 on the clean Ge4Ru4 data)
    @test maximum(opnorm(Ud[:, :, k]'Ud[:, :, k] - I) for k in 1:sc.nk_ibz) < 1.0e-5

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    ws2 = Wannier.SymmetricTransportWorkspace(Ef, frozen, sc)
    fg! = function (F, G, x)
        Wannier.assemble_gauge!(ws2.U_ibz, x, sb)
        Ω = Wannier._fg_transport_core!(F, G === nothing ? nothing : ws2.G_ibz, mmn_i.M, sc, ws2)
        G === nothing || Wannier.pullback_gradient!(G, ws2.G_ibz, x, sb)
        return Ω
    end
    g = zero(x0)
    Ω = fg!(1.0, g, x0)
    @test isfinite(Ω)

    # directional FD of the Schur objective (exact chain; tolerance set by
    # FD truncation at this dataset's Im-log conditioning)
    for _ in 1:2
        dx = randn(length(x0))
        dx ./= norm(dx)
        ε = 1.0e-4
        fd = (fg!(1.0, nothing, x0 .+ ε .* dx) - fg!(1.0, nothing, x0 .- ε .* dx)) / (2ε)
        an = real(sum(conj.(g) .* dx))
        @test isapprox(fd, an; rtol = 1.0e-3)
    end
end

@testitem "energy masking and breaking force" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")

    sc0 = Wannier.SymmetryConstraint(ks0, isym, centers)
    sc = Wannier.SymmetryConstraint(ks0, isym, centers; eig_ibz = Ei)
    # energy-based masking can only extend the d-norm-based mask, and it
    # masks whole degenerate clusters
    for iki in 1:sc.nk_ibz
        @test all(sc.covariant_bands[iki] .<= sc0.covariant_bands[iki])
        E = Ei[:, iki]
        for n in 1:(sc.nbands - 1)
            if E[n + 1] - E[n] <= 1.0e-4
                @test sc.covariant_bands[iki][n] == sc.covariant_bands[iki][n + 1]
            end
        end
    end

    # breaking force: finite, in (0, 1], and ~0 for the projected gradient
    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    ws2 = Wannier.SymmetricTransportWorkspace(Ef, frozen, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    U0 = Wannier.project_covariant(Wannier.lowdin_orthonormalize(Ai), sc)
    fb = Wannier.symmetry_breaking_force(U0, mmn_i.M, sc, ws2)
    @test 0 <= fb <= 1
end

@testitem "clean_littlegroup_reps!" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    raw = [copy(Matrix(r.d)) for r in isym.littlegroup_reps]
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    atol_deg, atol_unit = 1.0e-4, 0.01
    clean_littlegroup_reps!(
        isym.littlegroup_reps, Ei;
        atol_degeneracy = atol_deg, atol_unitary = atol_unit,
    )

    # per rep: exactly block-diagonal over energy multiplets; each block is
    # either exactly unitarized (was near-unitary up to data noise) or a
    # genuinely truncated multiplet kept bit-identical to the raw data
    n_unitarized, n_kept = let n_unitarized = 0, n_kept = 0
        for (ir, rep) in enumerate(isym.littlegroup_reps)
            E = Ei[:, rep.ik_ibz]
            nb = size(rep.d, 1)
            cluster = zeros(Int, nb)
            lo = 1
            for n in 1:nb
                cluster[n] = lo
                (n == nb || E[n + 1] - E[n] > atol_deg) && (lo = n + 1)
            end
            @test all(
                rep.d[m, n] == 0 for m in 1:nb, n in 1:nb if cluster[m] != cluster[n]
            )
            for blk in (findall(==(c), cluster) for c in unique(cluster))
                B = Matrix(rep.d[blk, blk])
                if opnorm(B' * B - I) <= 1.0e-12
                    n_unitarized += 1
                else
                    @test B == raw[ir][blk, blk]         # left untouched
                    @test opnorm(B' * B - I) > atol_unit # a true contraction
                    n_kept += 1
                end
            end
        end
        (n_unitarized, n_kept)
    end
    @test n_unitarized > 0
    @test n_kept > 0   # Si2_hse's window truncates some multiplets

    # the symmetry-broken-band masking is unaffected by the cleaning, and
    # still detects Si2_hse's truncated multiplets
    centers = [p.center for p in nnkp["projections"]]
    isym_raw = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym_raw.littlegroup_reps)
    sc_raw = Wannier.SymmetryConstraint(kstencil, isym_raw, centers)
    sc = Wannier.SymmetryConstraint(kstencil, isym, centers)
    @test sc.covariant_bands == sc_raw.covariant_bands
    @test any(iki -> !all(sc.covariant_bands[iki]), 1:sc.nk_ibz)

    # cleaning lowers the data-noise floors by ~100x (Si2_hse raw d noise is
    # ~1e-5): projector idempotence and the covariance residual of the
    # projected .iamn drop well below the raw floors
    X = randn(ComplexF64, sc.nbands, sc.nwann, sc.nk_ibz)
    PX_raw = Wannier.project_covariant(X, sc_raw)
    PX = Wannier.project_covariant(X, sc)
    @test norm(Wannier.project_covariant(PX_raw, sc_raw) - PX_raw) > 1.0e-6
    @test norm(Wannier.project_covariant(PX, sc) - PX) < 1.0e-6

    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    @test Wannier.covariance_residual(Wannier.project_covariant(Ai, sc_raw), sc_raw) > 1.0e-6
    @test Wannier.covariance_residual(Wannier.project_covariant(Ai, sc), sc) < 1.0e-7
end

@testitem "localize on SymmetricModel" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(ks0, isym, centers)
    ks = Wannier.globalize_bvector_ordering(ks0)

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Mf = Wannier.reconstruct_overlaps(mmn_i.M, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    atom_positions = [p.second for p in win["atoms_frac"]]
    atom_labels = map(x -> string(x.first), win["atoms_frac"])
    model = Wannier.Model(
        win["unit_cell_cart"], atom_positions, atom_labels,
        ks, Mf, Wannier.reconstruct_gauges(Wannier.project_covariant(Ai, sc), sc), Ef, frozen,
    )

    sm = SymmetricModel(model, sc, mmn_i.M)
    @test Wannier.n_bands(sm) == Wannier.n_bands(model)
    @test Wannier.n_wannier(sm) == Wannier.n_wannier(model)
    @test Wannier.n_kpoints(sm) == sc.nk_fbz
    @test Wannier.n_kpoints_ibz(sm) == sc.nk_ibz
    @test default_layout(Variance(), sm) == SymmetricXYLayout()
    @test SymmetricXYLayout().path == :transport

    # framework path (Problem + solve! under the hood), a few LBFGS iterations
    niter = 5
    U_fbz, U_ibz = localize(sm; max_iter = niter)
    # the optimized gauge is covariant to the data floor and consistently reconstructed
    @test Wannier.covariance_residual(U_ibz, sc) < 1.0e-3
    @test U_fbz ≈ Wannier.reconstruct_gauges(U_ibz, sc)

    # hard-coded 5-iteration regression anchor (SymmetricXYLayout, transport path)
    Ω = Wannier.spread(model.kstencil, model.overlaps, U_fbz).Ω
    @test isapprox(Ω, 22.43757274472156; atol = 1.0e-7)

    # SAWF improves on covariance-projected input and stays a semi-unitary gauge
    @test Ω < Wannier.spread(model.kstencil, model.overlaps, model.gauges).Ω
    @test maximum(
        opnorm(U_ibz[:, :, k]' * U_ibz[:, :, k] - I) for k in 1:sc.nk_ibz
    ) < 1.0e-4

    # Schur layout through the explicit Problem + solve! route
    prob = Problem(Variance(), sm, SchurLayout())
    Us_fbz, Us_ibz = solve!(prob, OptimLBFGS(; max_iter = niter))
    @test Wannier.covariance_residual(Us_ibz, sc) < 1.0e-3
    # hard-coded 5-iteration regression anchor (SchurLayout)
    Ωs = Wannier.spread(model.kstencil, model.overlaps, Us_fbz).Ω
    @test isapprox(Ωs, 22.07929818239718; atol = 1.0e-7)

    # full-mesh path through the layout: same variables, same optimum as transport
    U1_fbz, _ = localize(Variance(), sm, SymmetricXYLayout(:fullmesh); max_iter = niter)
    Ω1 = Wannier.spread(model.kstencil, model.overlaps, U1_fbz).Ω
    @test isapprox(Ω1, 22.43757275126563; atol = 1.0e-7)
    @test isapprox(Ω1, Ω; atol = 1.0e-5)
end

@testitem "CenteredVariance on SymmetricModel" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.SymmetryConstraint(ks0, isym, centers)
    ks = Wannier.globalize_bvector_ordering(ks0)

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Mf = Wannier.reconstruct_overlaps(mmn_i.M, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    atom_positions = [p.second for p in win["atoms_frac"]]
    atom_labels = map(x -> string(x.first), win["atoms_frac"])
    model = Wannier.Model(
        win["unit_cell_cart"], atom_positions, atom_labels,
        ks, Mf, Wannier.reconstruct_gauges(Wannier.project_covariant(Ai, sc), sc), Ef, frozen,
    )
    sm = SymmetricModel(model, sc, mmn_i.M)

    # target centers: those of the starting covariant gauge, shifted so the
    # penalty gradient is nonzero
    r0 = Wannier.spread(model.kstencil, model.overlaps, model.gauges).r
    r0 = [r + Wannier.Vec3(0.05, -0.03, 0.02) for r in r0]
    λ = 1.0
    obj = CenteredVariance(r0, λ)

    prob2 = Problem(obj, sm)                    # default: SymmetricXYLayout (transport path)
    prob1 = Problem(obj, sm, SymmetricXYLayout(:fullmesh))
    fg2 = Wannier._optimizer_callback(prob2)
    fg1 = Wannier._optimizer_callback(prob1)
    x = Wannier.initial_parameters(prob2.layout, sm)
    g1, g2 = zero(x), zero(x)
    Ω1 = fg1(1.0, g1, x)
    Ω2 = fg2(1.0, g2, x)

    # the full-mesh path equals the full-mesh penalized spread of the expanded covariant
    # gauge exactly; the transport path agrees to the isym data noise (as for Variance)
    Uf, _ = Wannier.finalize_result(prob2.layout, x, sm)
    Ωt_ref = Wannier.omega_center(model.kstencil, model.overlaps, Uf; r0, λ).Ωt
    @test isapprox(Ω1, Ωt_ref; atol = 1.0e-10)
    @test isapprox(Ω2, Ω1; atol = 1.0e-6)
    @test norm(g1 - g2) / norm(g1) < 1.0e-4

    # directional FD of the penalized transport-path objective (rtol limited by FD
    # truncation at this dataset's near-branch-point Im-log diagonals, which
    # the center penalty amplifies; cf. the 1e-3 of the Schur FD test)
    for _ in 1:2
        dx = randn(ComplexF64, size(x))
        dx ./= norm(dx)
        ε = 1.0e-4
        fd = (fg2(1.0, nothing, x .+ ε .* dx) - fg2(1.0, nothing, x .- ε .* dx)) / (2ε)
        an = real(sum(conj.(g2) .* dx))
        @test isapprox(fd, an; rtol = 1.0e-3)
    end

    # a few localize iterations stay covariant and beat the starting value
    Uc_fbz, Uc_ibz = localize(obj, sm; max_iter = 3)
    @test Wannier.covariance_residual(Uc_ibz, sc) < 1.0e-3
    Ωc = Wannier.omega_center(model.kstencil, model.overlaps, Uc_fbz; r0, λ).Ωt
    @test Ωc < Ωt_ref
end
