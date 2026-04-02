@testmodule CooptCenterEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier.Datasets
    export model, f, g!, λs, p

    # model_up = read_w90(dataset"Fe_collinear/Fe_up")
    # model_dn = read_w90(dataset"Fe_collinear/Fe_dn")
    # Mupdn = read_amn(joinpath(dataset"Fe_collinear/Fe_updn.mud"))
    model_up = read_w90(joinpath(@__DIR__, "../../fixtures/iron_dftk", "up"))
    model_dn = read_w90(joinpath(@__DIR__, "../../fixtures/iron_dftk", "dn"))
    Mupdn = read_amn(joinpath(@__DIR__, "../../fixtures/iron_dftk", "iron.mud"))
    model = Wannier.MagModel(model_up, model_dn, Mupdn)
    # if λs=0, equivalent to two independent Wannierizations of up and down
    # λs = 0
    λs = 1.0
    r₀ = [Wannier.Vec3(zeros(3)) for i in 1:n_wannier(model.up)]
    λc = 10.0
    p = CenterSpreadPenalty(r₀, λc)
    f, g! = Wannier.get_fg!_disentangle(p, model, λs)
end

@testitem "coopt center spread" setup=[CooptCenterEnv] begin
    Ω = Wannier.omega(p, model, λs)
    @test isapprox(Ω.up.Ω, 5.926927084598657; atol=1e-10)
    @test isapprox(Ω.up.Ωt, 193.35456420417805; atol=1e-10)
    @test isapprox(Ω.dn.Ω, 5.857837450506195; atol=1e-10)
    @test isapprox(Ω.dn.Ωt, 192.08254963648324; atol=1e-10)
    @test isapprox(Ω.Ωupdn, 0.02594482662336972; atol=1e-10)
    @test isapprox(Ω.Ωt, 385.4630586672846; atol=1e-10)

    M = [
        0.9959313044867574     6.25429314670281e-6    2.128513404763234e-7   1.330679062695095e-5   3.3410941167332526e-6  1.9116231862363496e-7
        1.0403595099641532e-5  0.9953197498405131     9.416016685007132e-7   8.977186962336146e-7   5.95480315312971e-6    3.475651745203065e-8
        2.7918925472446733e-7  1.3223079816980856e-7  0.9976876577962708     6.876745482878178e-7   4.1822715623681855e-6  3.2003662646437373e-7
        7.939195345282984e-6   3.127255097243016e-6   6.28031800550074e-7    0.9937844582468215     3.7913597967603723e-6  1.2711124421676069e-6
        5.605735018341171e-7   2.846907338937114e-6   2.8897771598613697e-8  4.420780051550379e-7   0.9933458780019229     3.908797712134863e-6
        1.7307795290532862e-6  1.204638950893186e-7   1.9677651181284534e-7  1.0247459022465658e-7  4.092959231915451e-6   0.9979861250043451
    ]
    @test isapprox(Ω.M, M; atol=1e-10)
end

@testitem "coopt center overlap gradient" setup=[CooptCenterEnv] begin
    using NLSolversBase

    function fup(Uup)
        return Wannier.omega_updn(
            model, [Uup[:, :, ik] for ik in 1:size(Uup, 3)], model.dn.gauges
        )
    end
    function fdn(Udn)
        return Wannier.omega_updn(
            model, model.up.gauges, [Udn[:, :, ik] for ik in 1:size(Udn, 3)]
        )
    end

    # analytical gradient
    Gup, Gdn = Wannier.omega_updn_grad(model, model.up.gauges, model.dn.gauges)
    Gup *= λs
    Gdn *= λs

    # finite diff gradient
    u_up0 = stack(model.up.gauges)
    d = OnceDifferentiable(fup, u_up0)
    Gup_ref = NLSolversBase.gradient!(d, u_up0)
    u_dn0 = stack(model.dn.gauges)
    d = OnceDifferentiable(fdn, u_dn0)
    Gdn_ref = NLSolversBase.gradient!(d, u_dn0)

    # I am using a looser tolerance here
    @test isapprox(
        stack(Gup),
        Gup_ref;
        atol=1e-6,
    )
    @test isapprox(
        stack(Gdn),
        Gdn_ref;
        atol=1e-6,
    )
end

@testitem "coopt center spread gradient" setup=[CooptCenterEnv] begin
    using NLSolversBase

    nb, nw = size(model.up.gauges[1])
    n_inner = nb * nw + nw^2  # size of XY at each k-point

    Xup0, Yup0 = Wannier.U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = Wannier.U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
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
    Uup, Udn = Wannier.disentangle(p, model, λs; max_iter=1)

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
