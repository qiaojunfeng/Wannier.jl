@testitem "spread" begin
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")
    wout = read_wout(dataset"Si2_coarse/outputs/Si2.wout")

    Ω = spread(model)

    @test Ω.Ω ≈ wout["Ωtotal"]
    @test Ω.ΩI ≈ wout["ΩI"]
    @test Ω.ΩOD ≈ wout["ΩOD"]
    @test Ω.ΩD ≈ wout["ΩD"]
    @test Ω.Ω̃ ≈ wout["ΩD"] + wout["ΩOD"]

    @test isapprox(Ω.ω, wout["spreads"]; atol = 1.0e-8)
    @test all(isapprox.(Ω.r, wout["centers"]; atol = 1.0e-6))
end

@testitem "spread gradient" begin
    using NLSolversBase
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")
    fg! = Wannier._make_optim_fg!(Wannier.Problem(Wannier.Variance(), model, Wannier.UGauge()))

    U = copy(model.gauges)
    G = zero(U)
    fg!(nothing, G, U)

    # Use finite difference as reference
    Uinit = deepcopy(U)
    d = NLSolversBase.OnceDifferentiable(
        x -> fg!(1.0, nothing, x), Uinit, real(zero(eltype(Uinit)))
    )
    G_ref = NLSolversBase.gradient!(d, U)

    @test isapprox(G, G_ref; atol = 1.0e-7)
end

@testitem "center" begin
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")
    wout = read_wout(dataset"Si2_coarse/outputs/Si2.wout")

    r = center(model)
    @test all(isapprox.(r, wout["centers"]; atol = 1.0e-6))
end
