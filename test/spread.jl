@testitem "spread" begin
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")
    wout = read_wout(dataset"Si2_coarse/outputs/Si2.wout")

    Ω = spread(model)

    @test Ω.Ω ≈ wout["Ωtotal"]
    @test Ω.ΩI ≈ wout["ΩI"]
    @test Ω.ΩOD ≈ wout["ΩOD"]
    @test Ω.ΩD ≈ wout["ΩD"]
    @test Ω.Ωtilde ≈ wout["ΩD"] + wout["ΩOD"]

    @test isapprox(Ω.ω, wout["spreads"]; atol = 1.0e-8)
    @test all(isapprox.(Ω.r, wout["centers"]; atol = 1.0e-6))
end

@testitem "spread gradient" begin
    using NLSolversBase
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")
    fg! = Wannier._optimizer_callback(Wannier.Problem(Wannier.Variance(), model, Wannier.ULayout()))

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

@testitem "imaglog_guided" begin
    using LinearAlgebra
    # principal branch at zero guide
    for z in (1.0 + 0.5im, -1.0 + 1.0e-3im, cis(3.0))
        @test Wannier.imaglog_guided(z, 0.0) == Wannier.imaglog(z)
    end
    # picks the branch closest to -θ: continuous across the cut
    z = cis(3.14)   # imaglog ≈ +3.14
    @test Wannier.imaglog_guided(z, 3.15) ≈ 3.14 - 2π atol = 1.0e-12
    @test Wannier.imaglog_guided(z, -3.15) ≈ 3.14 atol = 1.0e-12
    # guided value differs from principal by an exact multiple of 2π
    for θ in (-7.0, 2.0, 9.9)
        d = Wannier.imaglog_guided(z, θ) - Wannier.imaglog(z)
        @test isapprox(rem(d, 2π, RoundNearest), 0; atol = 1.0e-12)
    end
end
