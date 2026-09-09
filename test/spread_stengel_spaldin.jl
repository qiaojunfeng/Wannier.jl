@testitem "StengelSpaldinPenalty ordering" begin
    using Wannier.Datasets
    model = load_dataset("Si2_valence")
    p = StengelSpaldinPenalty(model.kstencil)

    nbvecs = Wannier.n_bvectors(model.kstencil)
    # slot `ib` must hold the same b-vector at every kpoint
    for ik in 1:Wannier.n_kpoints(model.kstencil)
        bvecs = Wannier.get_bvectors(model.kstencil, ik)
        @test sort(p.nnord[ik]) == collect(1:nbvecs)
        @test all(
            isapprox(bvecs[p.nnord[ik][ib]], p.kstencil.bvectors[ib]; atol = 1.0e-6)
                for ib in 1:nbvecs
        )
    end

    # built from a different stencil than it is used with -> hard error, not
    # silently mixed b-vector conventions
    other = Wannier.reorder(model.kstencil)
    @test_throws ErrorException omega(p, other, model.overlaps, model.gauges)
end

@testitem "StengelSpaldinPenalty spread" begin
    using Wannier.Datasets
    model = load_dataset("Si2_valence")
    p = StengelSpaldinPenalty(model.kstencil)

    Ωmv = omega(model)
    Ωss = omega(p, model)
    @test Ωss.method === :StengelSpaldin
    @test Ωmv.method === :MarzariVanderbilt

    # SS changes only the diagonal part
    @test Ωss.ΩI ≈ Ωmv.ΩI
    @test Ωss.ΩOD ≈ Ωmv.ΩOD
    @test !isapprox(Ωss.ΩD, Ωmv.ΩD; atol = 1.0e-8)

    # ΩD is a variance over kpoints, hence non-negative
    @test Ωss.ΩD >= -1.0e-12
    # Ω = Σ_nb w_b (1 - |S_nb|²) collapses the three parts. Two independently
    # coded expressions for the same number, unlike Ω == ΩI + ΩOD + ΩD which is
    # how the `Spread` is constructed.
    @test Ωss.Ω ≈ sum(Ωss.ω)

    # the centres pin the sign, the 1/N_k and the w_b * b product; on this grid
    # both functionals localise about the same bond centres
    @test all(isapprox.(Ωss.r, Ωmv.r; atol = 1.0e-3))

    # the k-averaged overlaps obey S[n, -b] = conj(S[n, b]), which the gradient
    # derivation relies on
    cache = Wannier.Cache(model.kstencil, model.overlaps, model.gauges)
    Wannier.compute_MU_UtMU!(cache, model.kstencil, model.overlaps, model.gauges)
    S = Wannier.compute_ss_overlaps(p, cache.UtMU, Wannier.n_wannier(model))
    ibrev = [
        findfirst(b2 -> isapprox(b2, -b; atol = 1.0e-6), p.kstencil.bvectors)
            for b in p.kstencil.bvectors
    ]
    @test S[:, ibrev] ≈ conj(S)
end

@testitem "StengelSpaldinPenalty gradient" begin
    using NLSolversBase
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")
    p = StengelSpaldinPenalty(model.kstencil)
    fg! = Wannier.get_fg!_maxloc(p, model)

    nb = n_bands(model)
    nw = n_wannier(model)
    nk = n_kpoints(model)
    U = [model.gauges[ik][ib, ic] for ib in 1:nb, ic in 1:nw, ik in 1:nk]
    G = zero(U)
    fg!(nothing, G, U)

    # Use finite difference as reference
    Uinit = deepcopy(U)
    d = NLSolversBase.OnceDifferentiable(
        x -> fg!(1.0, nothing, x), Uinit, real(zero(eltype(Uinit)))
    )
    G_ref = NLSolversBase.gradient!(d, U)

    @test isapprox(G, G_ref; atol = 1.0e-7)
end

@testitem "max_localize functional keyword" begin
    using Wannier.Datasets
    model = load_dataset("Si2_valence")

    # the one-call entry point
    Umin = max_localize(model; functional = :StengelSpaldin, max_iter = 5)
    # same thing spelled out
    p = StengelSpaldinPenalty(model)
    @test Umin ≈ max_localize(p, model; max_iter = 5)

    # the default is unchanged
    @test max_localize(model; max_iter = 5) ≈
        max_localize(SpreadPenalty(), model; max_iter = 5)

    @test Wannier.spread_penalty(:MarzariVanderbilt, model) isa SpreadPenalty
    @test Wannier.spread_penalty(:StengelSpaldin, model) isa StengelSpaldinPenalty
    @test_throws ErrorException Wannier.spread_penalty(:nonsense, model)
end

@testitem "StengelSpaldinPenalty max_localize" begin
    using LinearAlgebra
    using Wannier.Datasets
    model = load_dataset("Si2_valence")
    # the penalty accepts a `Model` directly, so `max_localize` is a two-liner
    p = StengelSpaldinPenalty(model)

    Ωi = omega(p, model)
    Umin = max_localize(p, model; max_iter = 30)
    Ωf = omega(p, model, Umin)

    @test Ωf.Ω < Ωi.Ω
    # the gauge stays on the unitary manifold
    @test all(isapprox(U' * U, I; atol = 1.0e-10) for U in Umin)
    # ΩI is gauge invariant, so minimising must not move it
    @test Ωf.ΩI ≈ Ωi.ΩI
end
