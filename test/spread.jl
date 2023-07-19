@testitem "spread" begin
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/reference/Si2.chk")
    # wout = read_wout(dataset"Si2/reference/Si2.wout")
    wout = read_wout(expanduser("~/git/WannierDatasets/datasets/Si2/reference/Si2.wout"))

    Ω = omega(model)

    @test Ω.Ω ≈ 23.583436326
    @test Ω.ΩI ≈ 16.228844399
    @test Ω.ΩOD ≈ 7.093383270
    @test Ω.ΩD ≈ 0.261208658
    @test Ω.Ω̃ ≈ 7.354591928

    @test isapprox(Ω.ω, wout.spreads; atol=1e-5)
    @test isapprox(Ω.r, wout.centers; atol=1e-5)
end

@testitem "spread gradient" begin
    using NLSolversBase
    fg! = Wannier.get_fg!_maxloc(model)

    n_bands = size(model.U[1], 1)
    n_wann = size(model.U[1], 2)
    n_kpts = length(model.U)
    U = [model.U[ik][ib, ic] for ib in 1:n_bands, ic in 1:n_wann, ik in 1:n_kpts]
    G = zero(U)
    fg!(nothing, G, U)

    # Use finite difference as reference
    Uinit = deepcopy(U)
    d = NLSolversBase.OnceDifferentiable(
        x -> fg!(1.0, nothing, x), Uinit, real(zero(eltype(Uinit)))
    )
    G_ref = NLSolversBase.gradient!(d, U)

    @test isapprox(G, G_ref; atol=1e-7)
end

@testitem "center" begin
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/reference/Si2.chk")
    # wout = read_wout(dataset"Si2/reference/Si2.wout")
    wout = read_wout(expanduser("~/git/WannierDatasets/datasets/Si2/reference/Si2.wout"))

    r = center(model)
    @test isapprox(r, wout.centers; atol=1e-5)
end
