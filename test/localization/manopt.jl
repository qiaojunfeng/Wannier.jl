@testitem "ManoptLBFGS matches OptimLBFGS on Variance+UGauge" begin
    using Wannier.Datasets
    using Manopt, Manifolds

    model = load_dataset("Si2_valence")

    U_optim  = localize(Wannier.Variance(), model; g_tol = 1.0e-8, max_iter = 500)
    prob     = Wannier.Problem(Wannier.Variance(), model)
    U_manopt = Wannier.solve!(prob, ManoptLBFGS(; g_tol = 1.0e-8, max_iter = 500))

    Ω_optim  = spread(model, U_optim).Ω
    Ω_manopt = spread(model, U_manopt).Ω

    # both backends should converge to the same local minimum to high precision
    @test isapprox(Ω_optim, Ω_manopt; atol = 1.0e-6)
end

@testitem "ManoptLBFGS errors on unsupported variants" begin
    using Wannier.Datasets
    using Manopt, Manifolds

    # XYGauge path (entangled) is not yet ported
    model = load_dataset("Si2")
    prob = Wannier.Problem(Wannier.Variance(), model)
    @test_throws ErrorException Wannier.solve!(prob, ManoptLBFGS())
end
