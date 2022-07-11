using NLSolversBase

# A reusable fixture for a model
model = read_seedname("$FIXTURE_PATH/silicon")

f, g! = Wannier.get_fg!_maxloc(model)

@testset "spread" begin
    Ω = omega(model.bvectors, model.M, model.A)

    @test Ω.Ω ≈ 63.52324857972985
    @test Ω.ΩI ≈ 49.83393864002448
    @test Ω.ΩOD ≈ 13.671560897063946
    @test Ω.ΩD ≈ 0.0
    @test Ω.Ω̃ ≈ 13.689309939705375

    ω_ref = [7.279461 8.159524 8.164187 8.158470 7.279461 8.159511 8.164164 8.158471]'
    @test isapprox(Ω.ω, ω_ref; atol=1e-5)

    r_ref = [
        1.349402 1.348540 1.349402 1.348539 0.000000 0.000862 -0.000000 0.000863
        1.349402 1.348897 1.348903 1.349402 -0.000000 0.000505 0.000499 0.000000
        1.349402 1.349402 1.348683 1.348677 -0.000000 -0.000000 0.000719 0.000725
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
