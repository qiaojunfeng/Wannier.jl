using NLSolversBase

# A reusable fixture for a model
model = read_seedname(joinpath(FIXTURE_PATH, "silicon/silicon"))

f, g! = Wannier.get_fg!_maxloc(model)

@testset "spread" begin
    # should be roughly the same as test/fixtures/silicon/silicon.wout
    Ω = omega(model)

    @test Ω.Ω ≈ 18.19526254033958
    @test Ω.ΩI ≈ 12.274068863079536
    @test Ω.ΩOD ≈ 5.843656489756686
    @test Ω.ΩD ≈ 0.07753718750335814
    @test Ω.Ω̃ ≈ 5.921193677260044

    ω_ref = [1.763839 2.443091 2.447297 2.443499 1.763839 2.443618 2.448038 2.442041]'
    @test isapprox(Ω.ω, ω_ref; atol=1e-5)

    r_ref = [
        1.349402 1.348405 1.349402 1.348402 0.000000 0.001000 -0.000000 0.000997
        1.349402 1.348821 1.348821 1.349402 -0.000000 0.000586 0.000576 0.000000
        1.349402 1.349402 1.348567 1.348564 0.000000 -0.000000 0.000834 0.000839
    ]
    @test isapprox(Ω.r, r_ref; atol=1e-5)
end

@testset "spread gradient" begin
    G = zero(model.A)
    g!(G, model.A)

    # Use finite difference as reference
    Ainit = model.A
    d = NLSolversBase.OnceDifferentiable(f, Ainit, real(zero(eltype(Ainit))))
    G_ref = NLSolversBase.gradient!(d, model.A)

    @test isapprox(G, G_ref; atol=1e-7)
end

@testset "center" begin
    r = center(model)

    r_ref = [
        1.349402 1.348405 1.349402 1.348402 0.000000 0.001000 -0.000000 0.000997
        1.349402 1.348821 1.348821 1.349402 -0.000000 0.000586 0.000576 0.000000
        1.349402 1.349402 1.348567 1.348564 0.000000 -0.000000 0.000834 0.000839
    ]

    @test isapprox(r, r_ref; atol=1e-5)
end
