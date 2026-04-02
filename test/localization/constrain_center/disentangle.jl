@testmodule DisCenterEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier: Vec3
    using Wannier.Datasets
    export model, fg!, p

    # model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")

    # A reusable fixture for a model
    FIXTURE_PATH = joinpath(@__DIR__, "../../fixtures")
    model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))
    # Note I shift a little bit to avoid the center constraint being zero
    δ = 0.1
    r₀ = [
        [Vec3(1.34940, 1.34940, 1.34940) for i in 1:(n_wannier(model) / 2)]
        [Vec3(δ, δ, δ) for i in 1:(n_wannier(model) / 2)]
    ]
    λ = 10.0
    p = CenterSpreadPenalty(r₀, λ)
    fg! = Wannier.get_fg!_disentangle(p, model)
end

@testitem "constraint center disentangle spread gradient" setup=[DisCenterEnv] begin
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
    @test isapprox(G, G_ref; atol=1e-6)

    # Test 2nd iteration
    U1 = Wannier.disentangle(p, model; max_iter=1)
    X, Y = Wannier.U_to_X_Y(U1, model.frozen_bands)
    XY = Wannier.X_Y_to_XY(X, Y)

    fg!(nothing, G, XY)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)
    Wannier.zero_froz_grad!(G_ref, model.frozen_bands)
    @test isapprox(G, G_ref; atol=1e-6)
end

@testitem "constraint center disentangle" setup=[DisCenterEnv] begin
    using Wannier: Vec3

    Umin = Wannier.disentangle(p, model; max_iter=4)
    Ω = Wannier.omega(p, model.kstencil, model.overlaps, Umin)

    @test Ω.Ω ≈ Ω.ΩI + Ω.Ω̃
    @test Ω.Ω̃ ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 18.368092500164014; atol=1e-7)
    @test isapprox(Ω.ΩI, 11.767244242580883; atol=1e-7)
    @test isapprox(Ω.ΩOD, 6.479314551713672; atol=1e-7)
    @test isapprox(Ω.ΩD, 0.12153370586945841; atol=1e-7)
    @test isapprox(Ω.Ω̃, 6.60084825758313; atol=1e-7)

    @test isapprox(
        Ω.ω,
        [
            1.709716007908809,
            2.4661068964146278,
            2.4718662674188296,
            2.4714054479442975,
            1.7395685569804915,
            2.502800139071625,
            2.5105486039612495,
            2.4960805804640596,
        ];
        atol=1e-7,
    )
    @test isapprox(
        Ω.r,
        [
            [1.3664290815080815, 1.36638676095681, 1.3664008302974457],
            [1.3734774753294627, 1.37368488285987, 1.3663714193162906],
            [1.3664189673302354, 1.373752246189131, 1.373530735253152],
            [1.3735406090103803, 1.3663647113226876, 1.3736156191573443],
            [0.04529565812816302, 0.044922879610335434, 0.04516814752840035],
            [0.06708089867081801, 0.06662128599309228, 0.05039143477233052],
            [0.05056920730828514, 0.06683583881766089, 0.06717500634921911],
            [0.06687922585704462, 0.05001921517720845, 0.06672964951997101],
        ];
        atol=1e-7,
    )

    @test Ω.Ωt ≈ Ω.Ω + Ω.Ωc
    @test isapprox(Ω.Ωc, 0.28261141089279834; atol=1e-7)
    @test isapprox(Ω.Ωt, 18.650703911056812; atol=1e-7)
    @test isapprox(
        Ω.ωc,
        [
            0.008675678956152481,
            0.014575094273678646,
            0.014649695273013859,
            0.01456966644827729,
            0.09032586255515693,
            0.04658815525023828,
            0.04620745061233267,
            0.04701980752394818,
        ];
        atol=1e-7,
    )
    @test isapprox(
        Ω.ωt,
        [
            1.7183916868649614,
            2.480681990688306,
            2.4865159626918434,
            2.485975114392575,
            1.8298944195356484,
            2.5493882943218633,
            2.556756054573582,
            2.543100387988008,
        ];
        atol=1e-7,
    )
    @test Ω.ωt ≈ Ω.ω + Ω.ωc
end
