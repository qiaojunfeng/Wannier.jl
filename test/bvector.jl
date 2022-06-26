@testset "get_bvectors" begin
    win = read_win("$FIXTURE_PATH/silicon.win")
    _, kpb_k, kpb_b = read_mmn("$FIXTURE_PATH/silicon.mmn")

    kpoints = win["kpoints"]
    recip_lattice = get_recip_lattice(win["unit_cell"])

    bvectors = get_bvectors(kpoints, recip_lattice)

    ref_bvecs = zeros(Float64, 3, 8)
    ref_bvecs = [
        -0.291017 0.291017 -0.291017 0.291017 -0.291017 -0.291017 0.291017 0.291017
        -0.291017 -0.291017 0.291017 -0.291017 0.291017 -0.291017 0.291017 0.291017
        0.291017 -0.291017 -0.291017 0.291017 0.291017 -0.291017 0.291017 -0.291017
    ]

    ref_weights = zeros(Float64, 8)
    fill!(ref_weights, 1.4759541565587924)

    @test bvectors.recip_lattice ≈ recip_lattice
    @test bvectors.kpoints ≈ kpoints
    @test bvectors.weights ≈ ref_weights
    @test isapprox(bvectors.bvectors, ref_bvecs, atol = 1e-5)
    @test bvectors.kpb_k ≈ kpb_k
    @test bvectors.kpb_b ≈ kpb_b
end
