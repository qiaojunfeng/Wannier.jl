@testset "get_bvectors" begin
    win = read_win("$FIXTURE_PATH/silicon.win")
    _, kpb_k, kpb_b = read_mmn("$FIXTURE_PATH/silicon.mmn")

    kpoints = win["kpoints"]
    recip_lattice = get_recip_lattice(win["unit_cell"])

    bvectors = get_bvectors(kpoints, recip_lattice)

    @test bvectors.recip_lattice ≈ recip_lattice
    @test bvectors.kpoints ≈ kpoints
    # @test bvectors.bvectors ≈ ref_bvecs
    # @test bvectors.weights ≈ ref_weights
    @test bvectors.kpb_k ≈ kpb_k
    for ik = 1:size(kpb_k, 2)
        bvectors.kpb_b[:, :, ik] ≉ kpb_b[:, :, ik] && begin
            println(ik)
            println(bvectors.kpb_b[:, :, ik])
            println(kpb_b[:, :, ik])
            error("$ik")
        end
    end
end
