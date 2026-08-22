@testmodule MaxlocEnv begin
    using Wannier
    using Wannier.Datasets
    export model, fg!

    # no disentanglement
    model = read_w90_with_chk(dataset"Si2_valence_coarse/Si2", dataset"Si2_valence_coarse/outputs/Si2.chk")

    fg! = Wannier._make_fg!(Wannier.Problem(Wannier.Variance(), model, Wannier.ULayout()))
end

@testitem "maxloc spread gradient" setup = [MaxlocEnv] begin
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
    U1 = Wannier.localize(model; max_iter = 1)

    fg!(nothing, G, U1)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), U1, zero(eltype(real(U1))))
    G_ref = NLSolversBase.gradient!(d, U1)
    @test isapprox(G, G_ref; atol = 1.0e-6)
end

@testitem "maxloc valence" setup = [MaxlocEnv] begin
    using Wannier.Datasets
    # reset initial gauge
    model.gauges .= read_amn_ortho(dataset"Si2_valence_coarse/Si2.amn")

    Umin = localize(model)
    Ω = spread(model.kstencil, model.overlaps, Umin)

    @test isapprox(Ω.Ω, 4.086818459; atol = 1.0e-7)
    @test isapprox(Ω.ΩI, 3.706376532; atol = 1.0e-7)
    @test isapprox(Ω.ΩOD, 0.380441928; atol = 1.0e-7)
    @test isapprox(Ω.Ωtilde, 0.3804419269999997; atol = 1.0e-7)
end
