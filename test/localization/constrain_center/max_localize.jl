@testmodule MaxlocCenterEnv begin
    using Wannier
    using Wannier: Vec3
    using Wannier.Datasets
    export model, fg!, obj

    model = read_w90(dataset"Si2_valence_coarse/Si2")
    r0 = [Vec3(0.0, 0.0, 0.0) for i in 1:n_wannier(model)]
    λ = 10.0
    obj = Wannier.CenteredVariance(r0, λ)
    fg! = Wannier._make_fg!(Wannier.Problem(obj, model, Wannier.ULayout()))
end

@testitem "constraint center maxloc spread gradient" setup = [MaxlocCenterEnv] begin
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
    U1 = Wannier.localize(obj, model; max_iter = 1)

    fg!(nothing, G, U1)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), U1, zero(eltype(real(U1))))
    G_ref = NLSolversBase.gradient!(d, U1)
    @test isapprox(G, G_ref; atol = 1.0e-6)
end

@testitem "constraint center maxloc valence" setup = [MaxlocCenterEnv] begin
    Umin = Wannier.localize(obj, model; max_iter = 4)
    Ω = Wannier.omega_center(
        Wannier.spread(model.kstencil, model.overlaps, Umin);
        r0 = obj.r0,
        λ = obj.λ,
    )

    # NOTE: these fixtures intentionally track historical main-branch values
    # at a fixed low iteration count.
    @test Ω.Ω ≈ Ω.ΩI + Ω.Ω̃
    @test Ω.Ω̃ ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 30.11589846146567; atol = 1.0e-7)
    @test isapprox(Ω.ΩI, 3.706376531801815; atol = 1.0e-7)
    @test isapprox(Ω.ΩOD, 6.166406390333415; atol = 1.0e-7)
    @test isapprox(Ω.ΩD, 20.24311553933044; atol = 1.0e-7)

    @test isapprox(
        Ω.ω,
        [7.18477647718601, 7.5174589128747575, 8.485999979449879, 6.927663091955029];
        atol = 1.0e-7,
    )
    @test isapprox(
        Ω.r,
        [
            [0.09848736571647632, 0.10105426078948893, -0.3397567359100563],
            [-0.15843672629787864, -0.7074236628984307, 0.7030443899918356],
            [-0.12106696606828928, -0.10699898366786562, -0.0775247363763475],
            [0.007623200427533831, -0.16940228161351428, 0.5054673027513011],
        ];
        atol = 1.0e-7,
    )

    @test Ω.Ωt ≈ Ω.Ω + Ω.Ωc
    @test isapprox(Ω.Ωc, 14.715367316738073; atol = 1.0e-7)
    @test isapprox(Ω.Ωt, 44.83126577820374; atol = 1.0e-7)
    @test isapprox(
        Ω.ωc,
        [1.353463644257367, 10.198218493676137, 0.32116077529158654, 2.8425244035129826];
        atol = 1.0e-7,
    )
    @test isapprox(
        Ω.ωt,
        [8.538240121443376, 17.715677406550896, 8.807160754741465, 9.770187495468011];
        atol = 1.0e-7,
    )
    @test Ω.ωt ≈ Ω.ω + Ω.ωc
end
