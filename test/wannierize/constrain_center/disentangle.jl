using NLSolversBase


# A reusable fixture for a model
# no disentanglement
model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))
r₀ = [[Vec3(1.34940, 1.34940, 1.34940) for i = 1:model.n_wann/2]; [Vec3(0.0,0.0,0.0) for i = 1:model.n_wann/2]]
λ = 10.0
fg! = Wannier.get_fg!_center_disentangle(model, r₀, λ)

@testset "constraint center disentangle spread gradient" begin
    U0 = deepcopy(model.U)

    # analytical gradient
    X, Y = Wannier.U_to_X_Y(U0, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)
    G = similar(XY)
    fg!(nothing, G, XY)

    # finite diff gradient
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)

    # The gradient for frozen bands need to be set as 0 explicitly
    Wannier.zero_froz_grad!(G_ref, model.frozen_bands)
    @test isapprox(G, G_ref; atol=1e-6)

    # Test 2nd iteration
    U1 = Wannier.disentangle_center(model, r₀, λ; max_iter=1)
    X, Y = Wannier.U_to_X_Y(U1, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)

    fg!(nothing, G, XY)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)
    Wannier.zero_froz_grad!(G_ref, model.frozen_bands)
    @test isapprox(G, G_ref; atol=1e-6)
end

@testset "constraint center disentangle" begin
    Umin = Wannier.disentangle_center(model, r₀, λ; max_iter=4);
    Ω = Wannier.omega_center(model.bvectors, model.M, Umin; r₀, λ)

    display(Ω)
    @test Ω.Ω ≈ Ω.ΩI + Ω.Ω̃
    @test Ω.Ω̃ ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 18.13023065207267; atol=1e-7)
    @test isapprox(Ω.ΩI, 11.672695160371777; atol=1e-7)
    @test isapprox(Ω.ΩOD, 6.341223286965608; atol=1e-7)
    @test isapprox(Ω.ΩD, 0.11631220473528359; atol=1e-7)
    @test isapprox(Ω.Ω̃, 6.457535491700892; atol=1e-7)

    @test isapprox(
        Ω.ω,
        [
            1.6720494155597025,
            2.460415452388486,
            2.4669444682373154,
            2.4657363554531777,
            1.6720635891342623,
            2.4637620130565874,
            2.4718181805922725,
            2.457441177650863,
        ];
        atol=1e-7,
    )
    @test isapprox(
        Ω.r,
        [
            Vec3(1.3493218953976438, 1.3493613581763424,  1.3493454033039636),
            Vec3(1.3487523304656146, 1.3491260454857306,  1.3493415966363236),
            Vec3(1.349320637608541, 1.3491584356558521,  1.348968696047302),
            Vec3(1.348788920937121, 1.349359045780774,  1.3489543861927666),
            Vec3(8.250336965638709e-5, 3.423661589664105e-5,  5.730256938007102e-5),
            Vec3(0.0006279664190044885, 0.00023781998769869842,  5.6624178504504735e-5),
            Vec3(8.189201302412916e-5, 0.0002630773899438576,  0.0004277217819759191),
            Vec3(0.0006247774390737991, 3.728918078729418e-5,  0.0004318068660840961)
        ];
        atol=1e-7,
    )

    @test Ω.Ωt ≈ Ω.Ω + Ω.Ωc
    @test isapprox(Ω.Ωc, 2.6352789712419767e-5; atol=1e-7)
    @test isapprox(Ω.Ωt, 18.130257004862383; atol=1e-7)
    @test isapprox(
        Ω.ωc,
        [
            1.0574318662876213e-7,
            4.979378545481902e-6,
            2.5067482115453843e-6,
            5.736665343582181e-6,
            1.126253633027044e-7,
            4.541064675376399e-6,
            2.5886193767347236e-6,
            5.781945009767707e-6,
        ];
        atol=1e-7,
    )
    @test isapprox(
        Ω.ωt,
        [
            1.6720495213028892,
            2.4604204317670315,
            2.466946974985527,
            2.4657420921185214,
            1.6720637017596256,
            2.4637665541212628,
            2.471820769211649,
            2.457446959595873,
        ];
        atol=1e-7,
    )
    @test Ω.ωt ≈ Ω.ω + Ω.ωc
end
