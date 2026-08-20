@testitem "parallel_transport valence" begin
    using Wannier.Datasets
    # the conduction bands poles are chosen randomly, skip test on conduction.
    model = read_w90(dataset"Si2_valence_coarse/Si2")
    Umin, _ = parallel_transport(model)

    Uref = read_amn(dataset"Si2_valence_coarse/outputs/ptg.amn").A
    @test isapprox(Umin, Uref; atol = 1.0e-5)

    ϵ0, ϵ1 = Wannier.compute_error(model, Umin)
    ϵ0_ref = 0.8166518514231456
    ϵ1_ref = 0.8914034179176087
    @test isapprox(ϵ0, ϵ0_ref; atol = 1.0e-5)
    @test isapprox(ϵ1, ϵ1_ref; atol = 1.0e-5)
end

@testitem "localize(ParallelTransport()) matches parallel_transport" begin
    using Wannier.Datasets
    model = read_w90(dataset"Si2_valence_coarse/Si2")
    U1 = localize(ParallelTransport(), model)
    U2, _ = parallel_transport(model)
    @test U1 ≈ U2
end

@testitem "parallel_transport neg coord" begin
    using Wannier.Datasets
    # Test PTG with a kgrid that have negative coordinates, i.e., -0.25 instead of 0.75
    # Before commit a1b05ae, the `parallel_transport` function will fail with at `index_bvector`.
    model = read_w90(dataset"GaAs_coarse/GaAs")
    # only 4 valence bands as an isolated manifold
    model = Wannier.truncate(model, 1:4, 1:4)

    Umin, _ = parallel_transport(model)

    Uref = read_amn(dataset"GaAs_coarse/outputs/GaAs.val.ptg.amn").A
    @test isapprox(Umin, Uref; atol = 1.0e-5)

    ϵ0, ϵ1 = Wannier.compute_error(model, Umin)

    ϵ0_ref = 0.20036010968611806
    ϵ1_ref = 0.20103804200989328

    @test isapprox(ϵ0, ϵ0_ref; atol = 1.0e-5)
    @test isapprox(ϵ1, ϵ1_ref; atol = 1.0e-5)
end
