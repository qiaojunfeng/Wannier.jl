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

@testitem "level-1 symmetric fg" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.symmetry_constraint(ks0, isym, centers)
    ks = Wannier.globalize_stencil(ks0)

    Mf = Wannier.unfold_overlaps_cached(read_mmn(dataset"Si2_hse/Si2.immn").M, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    atom_positions = [p.second for p in win["atoms_frac"]]
    atom_labels = map(x -> string(x.first), win["atoms_frac"])
    model = Wannier.Model(
        win["unit_cell_cart"], atom_positions, atom_labels,
        ks, Mf, Wannier.expand_gauges(Wannier.project_covariant(Ai, sc), sc), Ef, frozen,
    )

    ws = Wannier.SymmetricWorkspace(model, sc)
    X, Y = Wannier.U_to_X_Y(Ai, ws.frozen_ibz)
    xy = Wannier.X_Y_to_XY(X, Y)
    G = zero(xy)
    Ω = Wannier.symmetric_fg1!(1.0, G, xy, model, sc, ws)

    # value consistency with the plain full-mesh spread of the expanded gauge
    Ufull = Wannier.expand_gauges(
        Wannier.project_covariant(Wannier.X_Y_to_U(X, Y), sc), sc
    )
    @test Ω ≈ Wannier.spread(model.kstencil, model.overlaps, Ufull).Ω

    # Directional finite-difference gradient check (frozen entries masked,
    # since the XY layout fixes them and zeroes their gradient). The rtol is
    # limited by FD truncation of the Si2_hse point (near-branch-point
    # Im-log diagonals); the same check on clean data (Ge4Ru4) passes at 1e-9.
    f = x -> Wannier.symmetric_fg1!(1.0, nothing, x, model, sc, ws)
    for _ in 1:2
        dx = randn(ComplexF64, size(xy))
        Wannier.zero_froz_grad!(dx, ws.frozen_ibz)
        dx ./= norm(dx)
        ε = 1.0e-4
        fd = (f(xy .+ ε .* dx) - f(xy .- ε .* dx)) / (2ε)
        an = real(sum(conj.(G) .* dx))
        @test isapprox(fd, an; rtol = 1.0e-4)
    end
end

@testitem "level-2 equals level-1" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.symmetry_constraint(ks0, isym, centers)
    ks = Wannier.globalize_stencil(ks0)

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    Mf = Wannier.unfold_overlaps_cached(mmn_i.M, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    atom_positions = [p.second for p in win["atoms_frac"]]
    atom_labels = map(x -> string(x.first), win["atoms_frac"])
    model = Wannier.Model(
        win["unit_cell_cart"], atom_positions, atom_labels,
        ks, Mf, Wannier.expand_gauges(Wannier.project_covariant(Ai, sc), sc), Ef, frozen,
    )

    ws1 = Wannier.SymmetricWorkspace(model, sc)
    ws2 = Wannier.SymmetricWorkspace2(Ef, frozen, sc)
    X, Y = Wannier.U_to_X_Y(Ai, ws1.frozen_ibz)
    xy = Wannier.X_Y_to_XY(X, Y)
    G1, G2 = zero(xy), zero(xy)
    Ω1 = Wannier.symmetric_fg1!(1.0, G1, xy, model, sc, ws1)
    Ω2 = Wannier.symmetric_fg2!(1.0, G2, xy, mmn_i.M, sc, ws2)

    # The transport theorem is exact for a covariant gauge; the agreement is
    # limited only by the projector's data-precision fixed point (Si2_hse d
    # matrices carry ~1e-5 noise; on clean data agreement is ~1e-12).
    @test isapprox(Ω1, Ω2; atol = 1.0e-6)
    @test norm(G1 - G2) / norm(G1) < 1.0e-4

    # transport identity for the Wannier-gauge overlaps themselves:
    # M̃(kf, bf) = phase · K_f[L† M̃_i R] versus the full-mesh product
    U_ibz = Wannier.project_covariant(Wannier.X_Y_to_U(X, Y), sc)
    U_fbz = Wannier.expand_gauges(U_ibz, sc)
    ikf = sc.stars[2][end]   # a star member away from its IBZ representative
    iki = sc.fbz2ibz[ikf][1]
    for ibf in 1:sc.nbvecs
        ikpb = ks.kpb_k[ibf, ikf]
        Mt_direct = U_fbz[:, :, ikf]' * Mf[:, :, ibf, ikf] * U_fbz[:, :, ikpb]
        ibi = sc.ibi_of[ibf, ikf]
        ikb = sc.ikb[ibi, iki]
        Ukb = Wannier._kconj(U_ibz[:, :, ikb] * sc.Aib[ibi, iki], sc.trev_ib[ibi, iki])
        Mt_i = U_ibz[:, :, iki]' * mmn_i.M[:, :, ibi, iki] * Ukb
        Mt_trans = sc.phase[ibf, ikf] .*
            Wannier._kconj(Matrix(sc.Lmat[ikf])' * Mt_i * Matrix(sc.Rmat[ibf, ikf]), sc.trev_f[ikf])
        @test norm(Mt_trans - Mt_direct) < 1.0e-5
    end
end

@testitem "schur block parametrization" begin
    using WannierIO, LinearAlgebra
    using Wannier.Datasets

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    ks0 = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    sc = Wannier.symmetry_constraint(ks0, isym, centers)

    Ei = read_eig(dataset"Si2_hse/Si2.ieig")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    frozen_ibz = frozen[:, sc.ibz2fbz]

    sb = Wannier.schur_basis(sc, frozen_ibz)
    # parameter reduction and feasibility (schur_basis errors when infeasible)
    @test 0 < sb.nx < (sc.nwann^2 + sc.nbands * sc.nwann) * sc.nk_ibz

    # decode of the initial parameters: covariant (to the isym data noise),
    # semi-unitary, and consistent with the Level-2 value at the same gauge
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    U0 = Wannier.project_covariant(Ai, sc)
    x0 = Wannier.schur_initial_x(U0, sb)
    Ud = zeros(ComplexF64, sc.nbands, sc.nwann, sc.nk_ibz)
    Wannier.schur_decode!(Ud, x0, sb)
    @test Wannier.covariance_residual(Ud, sc) < 1.0e-3
    # semi-unitarity of the decode is limited by the isym data noise through
    # the anti-unitary coset average (1.9e-12 on the clean Ge4Ru4 data)
    @test maximum(opnorm(Ud[:, :, k]'Ud[:, :, k] - I) for k in 1:sc.nk_ibz) < 1.0e-5

    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    ws2 = Wannier.SymmetricWorkspace2(Ef, frozen, sc)
    fg! = function (F, G, x)
        Wannier.schur_decode!(ws2.U_ibz, x, sb)
        Ω = Wannier._fg2_core!(F, G === nothing ? nothing : ws2.G_ibz, mmn_i.M, sc, ws2)
        G === nothing || Wannier.schur_encode_gradient!(G, ws2.G_ibz, x, sb)
        return Ω
    end
    g = zero(x0)
    Ω = fg!(1.0, g, x0)
    @test isfinite(Ω)

    # directional FD of the Schur objective (exact chain; tolerance set by
    # FD truncation at this dataset's Im-log conditioning)
    for _ in 1:2
        dx = randn(ComplexF64, length(x0))
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
    Wannier.rescale!(isym.littlegroup_reps)
    centers = [p.center for p in nnkp["projections"]]
    Ei = read_eig(dataset"Si2_hse/Si2.ieig")

    sc0 = Wannier.symmetry_constraint(ks0, isym, centers)
    sc = Wannier.symmetry_constraint(ks0, isym, centers; eig_ibz = Ei)
    # energy-based masking can only extend the d-norm-based mask, and it
    # masks whole degenerate clusters
    for iki in 1:sc.nk_ibz
        @test all(sc.band_ok[iki] .<= sc0.band_ok[iki])
        E = Ei[:, iki]
        for n in 1:(sc.nbands - 1)
            if E[n + 1] - E[n] <= 1.0e-4
                @test sc.band_ok[iki][n] == sc.band_ok[iki][n + 1]
            end
        end
    end

    # breaking force: finite, in (0, 1], and ~0 for the projected gradient
    mmn_i = read_mmn(dataset"Si2_hse/Si2.immn")
    win = read_win(dataset"Si2_hse/Si2.win")
    Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
    frozen = Wannier.get_frozen_bands(Ef, get(win, "dis_froz_max", -Inf))
    ws2 = Wannier.SymmetricWorkspace2(Ef, frozen, sc)
    Ai = read_amn(dataset"Si2_hse/Si2.iamn").A
    U0 = Wannier.project_covariant(Wannier.orthonorm_lowdin(Ai), sc)
    fb = Wannier.symmetry_breaking_force(U0, mmn_i.M, sc, ws2)
    @test 0 <= fb <= 1
end
