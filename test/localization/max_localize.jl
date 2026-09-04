@testmodule MaxlocEnv begin
    using Wannier
    using Wannier.Datasets
    export model, fg!

    # no disentanglement
    model = read_w90_with_chk(dataset"Si2_valence_coarse/Si2", dataset"Si2_valence_coarse/outputs/Si2.chk")

    fg! = Wannier._optimizer_callback(Wannier.Problem(Wannier.Variance(), model, Wannier.ULayout()))
end

@testitem "maxloc spread gradient" setup = [MaxlocEnv] begin
    using NLSolversBase

    U0 = copy(model.gauges)
    # analytical gradient
    G = similar(U0)
    fg!(nothing, G, U0)

    # finite diff gradient
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), U0, zero(eltype(real(U0))))
    G_ref = NLSolversBase.gradient!(d, U0)

    # I am using a looser tolerance here
    @test isapprox(G, G_ref; atol = 1.0e-6)

    # Test 2nd iteration
    U1 = Wannier.localize(model; max_iter = 1)

    fg!(nothing, G, U1)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), U1, zero(eltype(real(U1))))
    G_ref = NLSolversBase.gradient!(d, U1)
    @test isapprox(G, G_ref; atol = 1.0e-6)
end

@testitem "maxloc valence" setup = [MaxlocEnv] begin
    using Wannier.Datasets
    # reset initial gauge
    model.gauges .= read_amn_ortho(dataset"Si2_valence_coarse/Si2.amn")

    Umin = localize(model)
    Ω = spread(model.kstencil, model.overlaps, Umin)

    @test isapprox(Ω.Ω, 4.086818459; atol = 1.0e-7)
    @test isapprox(Ω.ΩI, 3.706376532; atol = 1.0e-7)
    @test isapprox(Ω.ΩOD, 0.380441928; atol = 1.0e-7)
    @test isapprox(Ω.Ωtilde, 0.3804419269999997; atol = 1.0e-7)
end

@testitem "localization result" setup = [MaxlocEnv] begin
    problem = Wannier.Problem(Wannier.Variance(), model)
    result = solve(
        problem,
        OptimLBFGS(; max_iter = 2, store_trace = true, show_every = 0);
        warmup = true,
    )

    @test result isa LocalizationResult
    @test size(result.solution) == size(model.gauges)
    @test result.objective_value isa Float64
    @test result.gradient_norm >= 0
    @test result.iterations <= 2
    @test result.termination_reason in (
        :gradient_tolerance,
        :objective_tolerance,
        :parameter_tolerance,
        :iteration_limit,
        :line_search_failure,
        :unknown,
    )
    @test result.elapsed_seconds >= 0
    @test !isempty(result.trace)
    @test all(entry -> entry isa LocalizationTraceEntry, result.trace)
    @test localize(Wannier.Variance(), model; max_iter = 2) ≈ result.solution
    @test_throws ArgumentError OptimLBFGS(; show_every = -1)

    rendered = sprint(show, MIME"text/plain"(), result)
    @test occursin("termination_reason", rendered)
end
