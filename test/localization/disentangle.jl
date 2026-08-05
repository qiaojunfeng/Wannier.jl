@testmodule DisentangleEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier.Datasets
    export model, fg!

    model = read_w90(dataset"Si2_coarse/Si2")
    fg! = Wannier._make_fg!(Wannier.Problem(Wannier.Variance(), model, Wannier.XYLayout()))
end

@testitem "U_to_X_Y X_Y_to_U" setup = [DisentangleEnv] begin
    X, Y = Wannier.U_to_X_Y(model.gauges, model.frozen_bands)
    U1 = Wannier.X_Y_to_U(X, Y)
    # U1 != model.U since some states are frozen
    X1, Y1 = Wannier.U_to_X_Y(U1, model.frozen_bands)
    # X1 != X, Y1 != Y, since the X, Y gauge are arbitrary due to SVD
    # However the U1 should = U2
    U2 = Wannier.X_Y_to_U(X1, Y1)
    @test isapprox(U1, U2; atol = 1.0e-6)
end

@testitem "XY_to_X_Y X_Y_to_XY" setup = [DisentangleEnv] begin
    X, Y = Wannier.U_to_X_Y(model.gauges, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)
    X1, Y1 = Wannier.XY_to_X_Y(XY, n_bands(model), n_wannier(model))
    @test isapprox(X, X1; atol = 1.0e-6)
    @test isapprox(Y, Y1; atol = 1.0e-6)
end

@testitem "disentangle spread gradient" setup = [DisentangleEnv] begin
    using NLSolversBase

    U0 = deepcopy(model.gauges)

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
    @test isapprox(G, G_ref; atol = 1.0e-6)

    # Test 2nd iteration
    U1 = Wannier.localize(model; max_iter = 1)
    X, Y = Wannier.U_to_X_Y(U1, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)

    fg!(nothing, G, XY)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)
    Wannier.zero_froz_grad!(G_ref, model.frozen_bands)
    @test isapprox(G, G_ref; atol = 1.0e-6)
end

@testitem "disentangle" setup = [DisentangleEnv] begin
    Umin = Wannier.localize(model; max_iter = 4)
    Ω = Wannier.spread(model.kstencil, model.overlaps, Umin)

    # display(Ω)
    @test Ω.Ω ≈ Ω.ΩI + Ω.Ω̃
    @test Ω.Ω̃ ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 12.362335109447647; atol = 1.0e-7)
    @test isapprox(Ω.ΩI, 7.212573765139664; atol = 1.0e-7)
    @test isapprox(Ω.ΩOD, 4.929849448152594; atol = 1.0e-7)
    @test isapprox(Ω.ΩD, 0.2199118961553884; atol = 1.0e-7)

    @test isapprox(
        Ω.ω,
        [
            1.2145730431631943,
            1.655531487632798,
            1.655531511083343,
            1.6555315112860653,
            1.2145730436877367,
            1.6555314861440165,
            1.6555315134652924,
            1.6555315129852053,
        ];
        atol = 1.0e-7,
    )
    @test isapprox(
        Ω.r,
        [
            [9.386119783477958e-9, -4.499475172965251e-8, 2.4981053542113326e-8],
            [-4.567961841769885e-8, 1.3926293592272931e-8, -4.076140836146991e-8],
            [-6.90776277038549e-8, -3.595155572896913e-8, -5.210872698786784e-8],
            [6.506912168541962e-8, 1.564160769284076e-7, 3.341713422557596e-8],
            [1.3576325402117597, 1.3576324575774867, 1.3576324919985168],
            [1.3576325378090364, 1.357632474151124, 1.357632542136112],
            [1.357632344302077, 1.3576324993020412, 1.3576325171624346],
            [1.3576324580526173, 1.3576326584763287, 1.3576324408622935],
        ];
        atol = 1.0e-7,
    )
end
