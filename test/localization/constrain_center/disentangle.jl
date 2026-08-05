@testmodule DisCenterEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier: Vec3
    using Wannier.Datasets
    export model, fg!, obj

    model = read_w90(dataset"Si2_coarse/Si2")
    # Note I shift a little bit to avoid the center constraint being zero
    δ = 0.1
    # In Cartesian coordinates
    a1, a2 = Ref(model.lattice) .* model.atom_positions
    r0 = [
        [Vec3(a1 .+ δ) for i in 1:(n_wannier(model) / 2)]
        [Vec3(a2) for i in 1:(n_wannier(model) / 2)]
    ]
    λ = 10.0
    obj = Wannier.CenteredVariance(r0, λ)
    fg! = Wannier._make_fg!(Wannier.Problem(obj, model, Wannier.XYLayout()))
end

@testitem "constraint center disentangle spread gradient" setup = [DisCenterEnv] begin
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
    U1 = Wannier.localize(obj, model; max_iter = 1)
    X, Y = Wannier.U_to_X_Y(U1, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)

    fg!(nothing, G, XY)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)
    Wannier.zero_froz_grad!(G_ref, model.frozen_bands)
    @test isapprox(G, G_ref; atol = 1.0e-6)
end

@testitem "constraint center disentangle" setup = [DisCenterEnv] begin
    using Wannier: Vec3

    Umin = Wannier.localize(obj, model; max_iter = 4)
    Ω = Wannier.omega_center(
        Wannier.spread(model.kstencil, model.overlaps, Umin);
        r0 = obj.r0,
        λ = obj.λ,
    )

    @test Ω.Ω ≈ Ω.ΩI + Ω.Ω̃
    @test Ω.Ω̃ ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 12.547748075975804; atol = 1.0e-7)
    @test isapprox(Ω.ΩI, 7.26151163267514; atol = 1.0e-7)
    @test isapprox(Ω.ΩOD, 5.053835771792866; atol = 1.0e-7)
    @test isapprox(Ω.ΩD, 0.2324006715077962; atol = 1.0e-7)

    @test isapprox(
        Ω.ω,
        [
            1.2250392970929476,
            1.692208243017581,
            1.6922082883000362,
            1.6922082466702766,
            1.2287573318601899,
            1.6724421742524544,
            1.6724422425099998,
            1.6724422522723206,
        ];
        atol = 1.0e-7,
    )
    @test isapprox(
        Ω.r,
        [
            [0.07726278880875197, 0.07726279743376369, 0.07726277449770608],
            [0.07554436001670126, 0.07554436540363324, 0.08754854184598646],
            [0.08754854392819636, 0.07554435975053093, 0.07554437488566383],
            [0.07554437022981468, 0.08754854963404013, 0.07554437643984441],
            [1.3796988502546963, 1.3796988790522562, 1.3796988735755522],
            [1.374075575531994, 1.3740755579424486, 1.37120627300088],
            [1.37120626146246, 1.3740755587136024, 1.3740755590792224],
            [1.3740755679738366, 1.3712062621460146, 1.3740755498012156],
        ];
        atol = 1.0e-7,
    )

    @test Ω.Ωt ≈ Ω.Ω + Ω.Ωc
    @test isapprox(Ω.Ωc, 0.09240287146622556; atol = 1.0e-7)
    @test isapprox(Ω.Ωt, 12.640150947442029; atol = 1.0e-7)
    @test isapprox(
        Ω.ωc,
        [
            0.01550942576835644,
            0.013511952006648821,
            0.013511946980524579,
            0.013511939673883746,
            0.014607737408166122,
            0.0072499620092900076,
            0.007249953719815391,
            0.007249953899540431,
        ];
        atol = 1.0e-7,
    )
    @test isapprox(
        Ω.ωt,
        [
            1.240548722861304,
            1.70572019502423,
            1.7057202352805607,
            1.7057201863441602,
            1.243365069268356,
            1.6796921362617445,
            1.6796921962298152,
            1.679692206171861,
        ];
        atol = 1.0e-7,
    )
    @test Ω.ωt ≈ Ω.ω + Ω.ωc
end
