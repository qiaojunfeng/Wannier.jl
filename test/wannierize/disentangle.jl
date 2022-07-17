using NLSolversBase

# A reusable fixture for a model
model = read_seedname(joinpath(FIXTURE_PATH, "silicon/silicon"))
f, g! = Wannier.get_fg!_disentangle(model)

@testset "A_to_X_Y X_Y_to_A" begin
    X, Y = Wannier.A_to_X_Y(model.A, model.frozen_bands)
    A1 = Wannier.X_Y_to_A(X, Y)
    # A1 != model.A since some states are frozen
    X1, Y1 = Wannier.A_to_X_Y(A1, model.frozen_bands)
    # X1 != X, Y1 != Y, since the X, Y gauge are arbitrary due to SVD
    # However the A1 should = A2
    A2 = Wannier.X_Y_to_A(X1, Y1)
    @test isapprox(A1, A2; atol=1e-6)
end

@testset "XY_to_X_Y X_Y_to_XY" begin
    X, Y = Wannier.A_to_X_Y(model.A, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)
    X1, Y1 = Wannier.XY_to_X_Y(XY, model.n_bands, model.n_wann)
    @test isapprox(X, X1; atol=1e-6)
    @test isapprox(Y, Y1; atol=1e-6)
end

"""Set grad of frozen bands to 0"""
function zero_froz_grad!(G::Matrix, model::Wannier.Model)
    GX, GY = Wannier.XY_to_X_Y(G, model.n_bands, model.n_wann)
    for ik in 1:size(G, 2)
        idx_f = model.frozen_bands[:, ik]
        n_froz = count(idx_f)
        GY[idx_f, :, ik] .= 0
        GY[:, 1:n_froz, ik] .= 0
    end
    G .= Wannier.X_Y_to_XY(GX, GY)
    return nothing
end

@testset "disentangle spread gradient" begin
    A0 = deepcopy(model.A)

    # analytical gradient
    X, Y = Wannier.A_to_X_Y(A0, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)
    G = similar(XY)
    g!(G, XY)

    # finite diff gradient
    d = OnceDifferentiable(f, XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)

    # The gradient for frozen bands need to be set as 0 explicitly
    zero_froz_grad!(G_ref, model)
    @test isapprox(G, G_ref; atol=1e-6)

    # Test 2nd iteration
    A1 = Wannier.disentangle(model; max_iter=1)
    X, Y = Wannier.A_to_X_Y(A1, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)

    g!(G, XY)
    d = OnceDifferentiable(f, XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)
    zero_froz_grad!(G_ref, model)
    @test isapprox(G, G_ref; atol=1e-6)
end

@testset "disentangle" begin
    Amin = Wannier.disentangle(model; max_iter=4)
    Ω = Wannier.omega(model.bvectors, model.M, Amin)

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
            1.349540244604139 1.3494591284925146 1.3494989831767787
            1.3483994786002111 1.3489677557627042 1.3492175157250734
            1.349148846101609 1.349041631018211 1.3487479885732723
            1.348490312952771 1.3492694250911146 1.348720943361702
            -0.00013820836974004684 -5.965595567273871e-5 -9.848597753393869e-5
            0.0009659756244854475 0.00035987176361155105 0.00018872798640355847
            0.00026663589947916155 0.0004041165538806068 0.0006513524118075078
            0.0009449317726510505 0.00011446347253599799 0.000661121718585481
        ]';
        atol=1e-7,
    )
end
