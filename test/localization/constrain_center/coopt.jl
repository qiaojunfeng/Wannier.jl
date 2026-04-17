@testmodule CooptCenterEnv begin
    using Wannier
    using Wannier.Datasets
    export model, f, g!, λs, obj

    model_up = read_w90(dataset"Fe_collinear_coarse/Fe_up")
    model_dn = read_w90(dataset"Fe_collinear_coarse/Fe_dn")
    Mupdn = read_amn(dataset"Fe_collinear_coarse/Fe_updn.mud").A
    model = Wannier.SpinModel(model_up, model_dn, Mupdn)
    λs = 1.0
    r₀ = [Wannier.Vec3(zeros(3)) for i in 1:n_wannier(model.up)]
    λc = 10.0
    obj = Wannier.CenteredCoOptVariance(r₀, λc, λs)
    fg_bundle = Wannier._make_optim_fg!(Wannier.Problem(obj, model))
    f, g! = fg_bundle
end

@testitem "coopt center spread" setup = [CooptCenterEnv] begin
    Ω_up = Wannier.spread(model.up)
    Ω_dn = Wannier.spread(model.dn)
    Ω_up_c = Wannier.omega_center(Ω_up; r₀ = obj.r0, λ = obj.λ)
    Ω_dn_c = Wannier.omega_center(Ω_dn; r₀ = obj.r0, λ = obj.λ)
    M_overlap = Wannier.overlap_updn(model)
    Ωupdn = Wannier.omega_updn(M_overlap)
    Ωt = Ω_up_c.Ωt + Ω_dn_c.Ωt + λs * Ωupdn

    @test isapprox(Ω_up.Ω, 5.962059896476422; atol = 1.0e-10)
    @test isapprox(Ω_up_c.Ωt, 7.581508110391737; atol = 1.0e-10)
    @test isapprox(Ω_dn.Ω, 6.361552430790301; atol = 1.0e-10)
    @test isapprox(Ω_dn_c.Ωt, 8.350994167955465; atol = 1.0e-10)
    @test isapprox(Ωupdn, 0.6214634458453414; atol = 1.0e-10)
    @test isapprox(Ωt, 16.553965724192544; atol = 1.0e-10)

    M = [
        0.99979691600997 1.0804427231765106e-20 6.709979456174779e-20 1.4172867197472119e-19 5.4562115775280394e-21 3.447623847355787e-22 3.679953459008748e-22 5.755890238894293e-22 6.411200210202371e-22
        2.0718940894181966e-20 0.78044159382794 0.0038141762881387445 0.007810835232019771 2.4199220522653525e-14 1.6699747842546642e-15 1.798502108153958e-15 3.6159201532174775e-11 5.037130854559965e-13
        8.167412067350954e-20 0.0003740146127939681 0.7845515562818077 0.0007498142653320732 2.3975988575521443e-13 6.8367652205297344e-15 8.084806895439038e-15 4.271521555890511e-10 3.633135931098627e-14
        2.2503786034881527e-19 0.008998430059967972 0.0019555181347221183 0.8140158111574739 4.1229951307567996e-17 1.5203647588094592e-16 1.456275022849239e-16 5.087478797265483e-13 3.6218213924869516e-13
        5.862479936630014e-21 2.5927331055586313e-12 8.499979391385873e-14 2.0126733063220623e-12 0.9999488253401912 2.5931543120781252e-24 1.2582785101267568e-23 3.919172407780589e-24 3.9700965960696657e-23
        4.263249118993395e-22 2.185145338738584e-15 5.475650500018344e-15 1.910726135685062e-15 3.1211905651685806e-22 0.9999443420611075 2.962119335348447e-23 9.355076486527723e-23 6.3688155469348456e-24
        2.2537602988733836e-22 2.196793623239565e-15 1.1475610598914752e-14 2.413207819924146e-15 4.547996184750753e-24 5.937108011375737e-24 0.9999443420597673 4.419035431651825e-25 6.209456477418765e-23
        4.5195847621531553e-23 1.1930744532518502e-12 3.1372910651089115e-14 9.301697510058695e-13 2.6089707342721644e-23 4.538841507408542e-24 1.1418418775586505e-24 0.9999488252974881 8.581365868804793e-24
        2.2401277849582583e-22 6.808579760244394e-15 9.892324401437166e-15 4.285861909579054e-15 1.443019694963365e-22 1.7005975766179996e-23 5.073760340439049e-24 2.0913947642096764e-23 0.9999443421189137
    ]
    @test isapprox(M_overlap, M; atol = 1.0e-10)
end

@testitem "coopt center overlap gradient" setup = [CooptCenterEnv] begin
    using NLSolversBase

    function fup(Uup)
        return Wannier.omega_updn(model, Uup, model.dn.gauges)
    end
    function fdn(Udn)
        return Wannier.omega_updn(model, model.up.gauges, Udn)
    end

    Gup, Gdn = Wannier.omega_updn_grad(model, model.up.gauges, model.dn.gauges)
    Gup *= λs
    Gdn *= λs

    u_up0 = copy(model.up.gauges)
    d = OnceDifferentiable(fup, u_up0)
    Gup_ref = NLSolversBase.gradient!(d, u_up0)
    u_dn0 = copy(model.dn.gauges)
    d = OnceDifferentiable(fdn, u_dn0)
    Gdn_ref = NLSolversBase.gradient!(d, u_dn0)

    @test isapprox(Gup, Gup_ref; atol = 1.0e-6)
    @test isapprox(Gdn, Gdn_ref; atol = 1.0e-6)
end

@testitem "coopt center spread gradient" setup = [CooptCenterEnv] begin
    using NLSolversBase

    nb, nw = size(model.up.gauges, 1), size(model.up.gauges, 2)
    n_inner = nb * nw + nw^2

    Xup0, Yup0 = Wannier.U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = Wannier.U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
    XY0 = vcat(Wannier.X_Y_to_XY(Xup0, Yup0), Wannier.X_Y_to_XY(Xdn0, Ydn0))

    G = similar(XY0)
    g!(G, XY0)

    d = OnceDifferentiable(f, XY0, zero(real(eltype(XY0))))
    G_ref = NLSolversBase.gradient!(d, XY0)
    Gu = @view G_ref[1:n_inner, :]
    Gd = @view G_ref[(n_inner + 1):end, :]
    Wannier.zero_froz_grad!(Gu, model.up.frozen_bands)
    Wannier.zero_froz_grad!(Gd, model.dn.frozen_bands)

    @test isapprox(G, G_ref; atol = 1.0e-6)

    # Test 2nd iteration
    Uup, Udn = Wannier.localize(obj, model; max_iter = 1)

    Xup0, Yup0 = Wannier.U_to_X_Y(Uup, model.up.frozen_bands)
    Xdn0, Ydn0 = Wannier.U_to_X_Y(Udn, model.dn.frozen_bands)
    XY0 = vcat(Wannier.X_Y_to_XY(Xup0, Yup0), Wannier.X_Y_to_XY(Xdn0, Ydn0))

    g!(G, XY0)
    d = OnceDifferentiable(f, XY0, zero(real(eltype(XY0))))
    G_ref = NLSolversBase.gradient!(d, XY0)
    # I need to use @view so that zero_froz_grad! can change it inplace
    Gu = @view G_ref[1:n_inner, :]
    Gd = @view G_ref[(n_inner + 1):end, :]
    # The gradient for frozen bands need to be set as 0 explicitly
    Wannier.zero_froz_grad!(Gu, model.up.frozen_bands)
    Wannier.zero_froz_grad!(Gd, model.dn.frozen_bands)

    @test isapprox(G, G_ref; atol = 1.0e-6)
end
