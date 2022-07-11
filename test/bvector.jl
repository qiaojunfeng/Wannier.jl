@testset "get_bvectors" begin
    win = read_win(joinpath(FIXTURE_PATH, "silicon.win"))
    _, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "silicon.mmn"))

    kpoints = win.kpoints
    recip_lattice = get_recip_lattice(win.unit_cell)

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
    @test begin
        # sometimes the bvectors are not ordered
        ret = true
        bvecs = Vector([bvectors.bvectors[:, i] for i in 1:size(bvectors.bvectors, 2)])
        for i in 1:size(ref_bvecs, 2)
            b = ref_bvecs[:, i]
            idx = findfirst(v -> isapprox(v, b; atol=1e-5), bvecs)
            if idx === nothing
                ret = false
                break
            end
            deleteat!(bvecs, idx)
        end
        length(bvecs) != 0 && (ret = false)
        ret
    end
    @test bvectors.kpb_k ≈ kpb_k
    @test bvectors.kpb_b ≈ kpb_b
end
