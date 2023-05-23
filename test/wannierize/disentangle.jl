using NLSolversBase

# A reusable fixture for a model
model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))
f, g! = Wannier.get_fg!_disentangle(model)

@testset "U_to_X_Y X_Y_to_U" begin
    X, Y = Wannier.U_to_X_Y(model.U, model.frozen_bands)
    U1 = Wannier.X_Y_to_U(X, Y)
    # U1 != model.U since some states are frozen
    X1, Y1 = Wannier.U_to_X_Y(U1, model.frozen_bands)
    # X1 != X, Y1 != Y, since the X, Y gauge are arbitrary due to SVD
    # However the U1 should = U2
    U2 = Wannier.X_Y_to_U(X1, Y1)
    @test isapprox(U1, U2; atol=1e-6)
end

@testset "XY_to_X_Y X_Y_to_XY" begin
    X, Y = Wannier.U_to_X_Y(model.U, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)
    X1, Y1 = Wannier.XY_to_X_Y(XY, model.n_bands, model.n_wann)
    @test isapprox(X, X1; atol=1e-6)
    @test isapprox(Y, Y1; atol=1e-6)
end

@testset "disentangle spread gradient" begin
    U0 = deepcopy(model.U)

    # analytical gradient
    X, Y = Wannier.U_to_X_Y(U0, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)
    G = similar(XY)
    g!(G, XY)

    # finite diff gradient
    d = OnceDifferentiable(f, XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)

    # The gradient for frozen bands need to be set as 0 explicitly
    Wannier.zero_froz_grad!(G_ref, model.frozen_bands)
    @test isapprox(G, G_ref; atol=1e-6)

    # Test 2nd iteration
    U1 = Wannier.disentangle(model; max_iter=1)
    X, Y = Wannier.U_to_X_Y(U1, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)

    g!(G, XY)
    d = OnceDifferentiable(f, XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)
    Wannier.zero_froz_grad!(G_ref, model.frozen_bands)
    @test isapprox(G, G_ref; atol=1e-6)
end

@testset "disentangle" begin
    Umin = Wannier.disentangle(model; max_iter=4)
    Ω = Wannier.omega(model.bvectors, model.M, Umin)

    # display(Ω)
    @test Ω.Ω ≈ Ω.ΩI + Ω.Ω̃
    @test Ω.Ω̃ ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 18.13021209358854; atol=1e-7)
    @test isapprox(Ω.ΩI, 11.672692248601113; atol=1e-7)
    @test isapprox(Ω.ΩOD, 6.341211889610597; atol=1e-7)
    @test isapprox(Ω.ΩD, 0.11630795537683003; atol=1e-7)
    @test isapprox(Ω.Ω̃, 6.457519844987427; atol=1e-7)

    @test isapprox(
        Ω.ω,
        [
            1.6720402167866713,
            2.4604144313827963,
            2.466947090855194,
            2.4657348106032355,
            1.672054382451226,
            2.463760873012476,
            2.471821111520987,
            2.457439176975955,
        ];
        atol=1e-7,
    )
    @test isapprox(
        Ω.r,
        [
            Vec3(1.349540244604139, 1.3494591284925146, 1.3494989831767787),
            Vec3(1.3483994786002111, 1.3489677557627042, 1.3492175157250734),
            Vec3(1.349148846101609, 1.349041631018211, 1.3487479885732723),
            Vec3(1.348490312952771, 1.3492694250911146, 1.348720943361702),
            Vec3(-0.00013820836974004684 , 5.965595567273871e-5, -9.848597753393869e-5),
            Vec3(0.0009659756244854475, 0.00035987176361155105, 0.00018872798640355847),
            Vec3(0.00026663589947916155, 0.0004041165538806068, 0.0006513524118075078),
            Vec3(0.0009449317726510505, 0.00011446347253599799, 0.000661121718585481)
        ];
        atol=1e-7,
    )
end
