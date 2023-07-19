using NLSolversBase
using LinearAlgebra

# A reusable fixture for a model
# no disentanglement
model_up = read_w90(joinpath(FIXTURE_PATH, "iron", "up"))
model_dn = read_w90(joinpath(FIXTURE_PATH, "iron", "dn"))
Mupdn = read_amn(joinpath(FIXTURE_PATH, "iron", "iron.mud"))
model = Wannier.MagModel(model_up, model_dn, Mupdn)
# if λ=0, equivalent to two independent Wannierizations of up and down
# λ = 0
λ = 1.0
f, g! = Wannier.get_fg!_disentangle(model, λ)

@testitem "coopt spread" begin
    Ω = Wannier.omega(model, λ)
    @test isapprox(Ω.up.Ω, 5.926926709743369; atol=1e-10)
    @test isapprox(Ω.dn.Ω, 5.857837078056738; atol=1e-10)
    @test isapprox(Ω.Ωupdn, 0.025659557291223933; atol=1e-10)
    @test isapprox(Ω.Ωt, 11.81042334509133; atol=1e-10)

    M = [
        0.9959665096311021 4.046700688326357e-5 2.185449243145392e-5 4.986028083817444e-5 3.934578618045495e-5 3.5746603114970393e-6
        5.048633998609008e-5 0.9953726121023302 2.6295300083381804e-5 6.910981438430957e-5 4.944071311188473e-5 1.1984326993369665e-5
        1.7943799040393165e-5 1.5478173490333126e-5 0.9976948923797521 2.558765638355381e-5 2.9520633577611506e-5 4.38144108797382e-6
        4.2778696242211456e-5 2.3135207831071605e-5 1.7669561968282892e-5 0.9939387319127542 8.221722397002391e-5 1.4724946905196543e-5
        3.3728247697111404e-5 5.4063099676496685e-5 2.2750384405991963e-5 4.7915516179350034e-5 0.993378311235527 2.6462003506350468e-5
        7.4369355742295675e-6 1.9465032795654377e-5 6.8005353881768465e-6 2.4244570561809903e-5 3.2040503339332296e-5 0.9979893854473111
    ]
    @test isapprox(Ω.M, M; atol=1e-10)
end

@testitem "coopt overlap gradient" begin
    n_bands, n_wann = size(model.up.U[1])
    n_kpts = length(model.up.U)

    function fup(Uup)
        return Wannier.omega_updn(
            model, [Uup[:, :, ik] for ik in 1:size(Uup, 3)], model.dn.U
        )
    end
    function fdn(Udn)
        return Wannier.omega_updn(
            model, model.up.U, [Udn[:, :, ik] for ik in 1:size(Udn, 3)]
        )
    end

    # analytical gradient
    Gup, Gdn = Wannier.omega_updn_grad(model, model.up.U, model.dn.U)
    Gup *= λ
    Gdn *= λ

    # finite diff gradient
    u_up0 = [
        model.up.U[ik][ib, ic] for ib in 1:size(model.up.U[1], 1),
        ic in 1:size(model.up.U[1], 2), ik in 1:length(model.up.U)
    ]
    d = OnceDifferentiable(fup, u_up0)
    Gup_ref = NLSolversBase.gradient!(d, u_up0)
    u_dn0 = [
        model.dn.U[ik][ib, ic] for ib in 1:size(model.dn.U[1], 1),
        ic in 1:size(model.dn.U[1], 2), ik in 1:length(model.dn.U)
    ]
    d = OnceDifferentiable(fdn, u_dn0)
    Gdn_ref = NLSolversBase.gradient!(d, u_dn0)

    # I am using a looser tolerance here
    @test isapprox(
        [
            Gup[ik][ib, ic] for ib in 1:size(Gup[1], 1), ic in 1:size(Gup[1], 2),
            ik in 1:length(Gup)
        ],
        Gup_ref;
        atol=1e-6,
    )
    @test isapprox(
        [
            Gdn[ik][ib, ic] for ib in 1:size(Gdn[1], 1), ic in 1:size(Gdn[1], 2),
            ik in 1:length(Gdn)
        ],
        Gdn_ref;
        atol=1e-6,
    )
end

@testitem "coopt spread gradient" begin
    n_bands, n_wann = size(model.up.U[1])
    n_kpts = length(model.up.U)
    n_inner = n_bands * n_wann + n_wann^2  # size of XY at each k-point

    Xup0, Yup0 = Wannier.U_to_X_Y(model.up.U, model.up.frozen_bands)
    Xdn0, Ydn0 = Wannier.U_to_X_Y(model.dn.U, model.dn.frozen_bands)
    # compact storage
    XYup0 = Wannier.X_Y_to_XY(Xup0, Yup0)
    XYdn0 = Wannier.X_Y_to_XY(Xdn0, Ydn0)
    XY0 = vcat(XYup0, XYdn0)

    # analytical gradient
    G = similar(XY0)
    g!(G, XY0)

    # finite diff gradient
    d = OnceDifferentiable(f, XY0, zero(real(eltype(XY0))))
    G_ref = NLSolversBase.gradient!(d, XY0)
    # I need to use @view so that zero_froz_grad! can change it inplace
    Gu = @view G_ref[1:n_inner, :]
    Gd = @view G_ref[(n_inner + 1):end, :]
    # The gradient for frozen bands need to be set as 0 explicitly
    Wannier.zero_froz_grad!(Gu, model.up.frozen_bands)
    Wannier.zero_froz_grad!(Gd, model.dn.frozen_bands)

    # I am using a looser tolerance here
    @test isapprox(G, G_ref; atol=1e-6)

    # Test 2nd iteration
    Uup, Udn = Wannier.disentangle(model, λ; max_iter=1)

    Xup0, Yup0 = Wannier.U_to_X_Y(Uup, model.up.frozen_bands)
    Xdn0, Ydn0 = Wannier.U_to_X_Y(Udn, model.dn.frozen_bands)
    # compact storage
    XYup0 = Wannier.X_Y_to_XY(Xup0, Yup0)
    XYdn0 = Wannier.X_Y_to_XY(Xdn0, Ydn0)
    XY0 = vcat(XYup0, XYdn0)

    g!(G, XY0)
    d = OnceDifferentiable(f, XY0, zero(real(eltype(XY0))))
    G_ref = NLSolversBase.gradient!(d, XY0)
    # I need to use @view so that zero_froz_grad! can change it inplace
    Gu = @view G_ref[1:n_inner, :]
    Gd = @view G_ref[(n_inner + 1):end, :]
    # The gradient for frozen bands need to be set as 0 explicitly
    Wannier.zero_froz_grad!(Gu, model.up.frozen_bands)
    Wannier.zero_froz_grad!(Gd, model.dn.frozen_bands)

    # I am using a looser tolerance here
    @test isapprox(G, G_ref; atol=1e-6)
end
