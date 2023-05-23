using Wannier: Vec3

@testset "get_bvectors" begin
    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))
    _, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "silicon/silicon.mmn"))

    kpoints = win.kpoints
    recip_lattice = get_recip_lattice(win.unit_cell_cart)

    bvectors = get_bvectors(kpoints, recip_lattice)

    ref_bvecs = [
        -0.291017 0.291017 -0.291017 0.291017 -0.291017 -0.291017 0.291017 0.291017
        -0.291017 -0.291017 0.291017 -0.291017 0.291017 -0.291017 0.291017 0.291017
        0.291017 -0.291017 -0.291017 0.291017 0.291017 -0.291017 0.291017 -0.291017
    ]
    ref_bvecs = map(i -> Vec3(ref_bvecs[:, i]), axes(ref_bvecs,2))

    ref_weights = zeros(Float64, 8)
    fill!(ref_weights, 1.4759541565587924)

    @test bvectors.recip_lattice ≈ recip_lattice
    @test bvectors.kpoints ≈ kpoints
    @test bvectors.weights ≈ ref_weights
    @test begin
        # sometimes the bvectors are not ordered
        ret = true
        bvecs = bvectors.bvectors
        for i in 1:length(ref_bvecs)
            b = ref_bvecs[i]
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

@testset "get_bvectors 2D" begin
    win = read_win(joinpath(FIXTURE_PATH, "graphene", "graphene.win"))
    nnkp = read_nnkp(joinpath(FIXTURE_PATH, "graphene", "graphene.nnkp"))

    kpoints = win.kpoints
    recip_lattice = get_recip_lattice(win.unit_cell_cart)

    bvectors = get_bvectors(kpoints, recip_lattice)

    ref_bvecs =
        [
            Vec3(0.000000, 0.070948, 0.000000),
            Vec3(0.061443, 0.035474, 0.000000),
            Vec3(-0.061443, -0.035474, 0.000000),
            Vec3(-0.061443, 0.035474, 0.000000),
            Vec3(0.000000, -0.070948, 0.000000),
            Vec3(0.061443, -0.035474, 0.000000),
            Vec3(0.000000, 0.000000, 0.209440),
            Vec3(0.000000, 0.000000, -0.209440),
        ]

    ref_weights = [
        66.220759
        66.220759
        66.220759
        66.220759
        66.220759
        66.220759
        11.398633
        11.398633
    ]

    @test bvectors.recip_lattice ≈ recip_lattice
    @test bvectors.kpoints ≈ kpoints
    @test bvectors.weights ≈ ref_weights
    @test begin
        # sometimes the bvectors are not ordered
        ret = true
        bvecs = [bvectors.bvectors[i] for i in 1:length(bvectors.bvectors)]
        for i in 1:length(ref_bvecs)
            b = ref_bvecs[i]
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
    @test bvectors.kpb_k ≈ nnkp.kpb_k
    @test bvectors.kpb_b ≈ nnkp.kpb_b
end

@testset "get_bvectors kmesh_tol" begin
    win = read_win(joinpath(FIXTURE_PATH, "kmesh_tol", "Se2Sn.win"))
    nnkp = read_nnkp(joinpath(FIXTURE_PATH, "kmesh_tol", "Se2Sn.nnkp"))

    kpoints = win.kpoints
    recip_lattice = get_recip_lattice(win.unit_cell_cart)

    bvectors = get_bvectors(kpoints, recip_lattice; win.kmesh_tol)

    ref_bvecs =
        [
            Vec3(0.000000, 0.000000, 0.180597),
            Vec3(0.000000, 0.000000, -0.180597),
            Vec3(0.188176, 0.000000, 0.000001),
            Vec3(-0.188176, 0.000000, -0.000001),
            Vec3(-0.094088, 0.162971, -0.000000),
            Vec3(0.094088, 0.162971, 0.000000),
            Vec3(0.094088, -0.162971, 0.000000),
            Vec3(-0.094088, -0.162971, -0.000000),
            Vec3(-0.188176, 0.000000, 0.180597),
            Vec3(0.188176, 0.000000, -0.180597)
        ]

    ref_weights = [
        15.330158
        15.330158
        9.413752
        9.413752
        9.412853
        9.412853
        9.412853
        9.412853
        0.000048
        0.000048
    ]

    @test bvectors.recip_lattice ≈ recip_lattice
    @test bvectors.kpoints ≈ kpoints
    @test all(isapprox.(bvectors.weights, ref_weights; atol=1e-5))
    @test begin
        # sometimes the bvectors are not ordered
        ret = true
        bvecs = [bvectors.bvectors[i] for i in 1:length(bvectors.bvectors)]
        for i in 1:length(ref_bvecs)
            b = ref_bvecs[i]
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
    @test bvectors.kpb_k ≈ nnkp.kpb_k
    @test bvectors.kpb_b ≈ nnkp.kpb_b
end
