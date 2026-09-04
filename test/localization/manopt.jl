@testitem "ManoptLBFGS matches OptimLBFGS on Variance+ULayout" begin
    using Wannier.Datasets
    using Manopt, Manifolds

    model = load_dataset("Si2_valence")

    U_optim = localize(Wannier.Variance(), model; g_tol = 1.0e-8, max_iter = 500)
    prob = Wannier.Problem(Wannier.Variance(), model)
    result = solve(
        prob,
        ManoptLBFGS(; g_tol = 1.0e-8, max_iter = 500, store_trace = true);
        warmup = true,
    )
    U_manopt = result.solution

    Ω_optim = spread(model, U_optim).Ω
    Ω_manopt = spread(model, U_manopt).Ω

    # both backends should converge to the same local minimum to high precision
    @test isapprox(Ω_optim, Ω_manopt; atol = 1.0e-6)
    @test result isa LocalizationResult
    @test result.converged
    @test result.termination_reason == :gradient_tolerance
    @test !isempty(result.trace)
end

@testitem "ManoptLBFGS errors on unsupported variants" begin
    using Wannier.Datasets
    using Manopt, Manifolds

    # XYLayout path (entangled) is not yet ported
    model = load_dataset("Si2")
    prob = Wannier.Problem(Wannier.Variance(), model)
    @test_throws ErrorException solve(prob, ManoptLBFGS())
end
