@testitem "generate_kspace_stencil" begin
    using Wannier.Datasets
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    mmn = read_mmn(dataset"Si2_valence/Si2_valence.mmn")
    kpb_k = mmn.kpb_k
    kpb_G = mmn.kpb_G

    recip_lattice = reciprocal_lattice(win["unit_cell_cart"])
    kstencil = generate_kspace_stencil(recip_lattice, win["mp_grid"], win["kpoints"])

    # copied from wout
    ref_bvectors = vec3.([
        [0.192835, 0.192835, -0.192835],
        [0.192835, -0.192835, 0.192835],
        [-0.192835, 0.192835, 0.192835],
        [0.192835, 0.192835, 0.192835],
        [-0.192835, -0.192835, 0.192835],
        [-0.192835, 0.192835, -0.192835],
        [0.192835, -0.192835, -0.192835],
        [-0.192835, -0.192835, -0.192835],
    ])
    ref_bweights = fill(3.361532, 8)
    ref_kstencil = Wannier.KspaceStencil(
        recip_lattice,
        vec3(win["mp_grid"]),
        win["kpoints"],
        ref_bvectors,
        ref_bweights,
        kpb_k,
        kpb_G,
    )
    @test isapprox(kstencil, ref_kstencil; atol = 1.0e-6)
end

@testitem "generate_kspace_stencil 2D" begin
    using Wannier.Datasets
    win = read_win(dataset"graphene_coarse/graphene.win")
    nnkp = read_nnkp_compute_bweights(dataset"graphene_coarse/outputs/graphene.nnkp")

    recip_lattice = reciprocal_lattice(win["unit_cell_cart"])
    kstencil = generate_kspace_stencil(recip_lattice, win["mp_grid"], win["kpoints"])

    ref_bvectors = vec3.([
        [0.0, 0.0, 0.628319],
        [0.0, 0.0, -0.628319],
        [0.793031, 0.457857, 0.0],
        [-0.793031, -0.457857, 0.0],
        [0.793031, -0.457857, 0.0],
        [-0.793031, 0.457857, 0.0],
        [0.0, -0.915713, 0.0],
        [0.0, 0.915713, 0.0],
    ])
    ref_bweights = [
        1.266515
        1.266515
        0.397521
        0.397521
        0.397521
        0.397521
        0.397521
        0.397521
    ]
    ref_kstencil = Wannier.KspaceStencil(
        recip_lattice,
        vec3(win["mp_grid"]),
        win["kpoints"],
        ref_bvectors,
        ref_bweights,
        nnkp.kpb_k,
        nnkp.kpb_G,
    )
    # wout does not have enough digits, use a bit larger atol
    @test isapprox(kstencil, ref_kstencil; atol = 1.0e-5)
end

@testitem "generate_kspace_stencil kmesh_tol" begin
    using Wannier.Datasets
    win = read_win(dataset"SnSe2/SnSe2.win")
    nnkp = read_nnkp_compute_bweights(dataset"SnSe2/outputs/SnSe2.nnkp")
    recip_lattice = reciprocal_lattice(win["unit_cell_cart"])
    kstencil = generate_kspace_stencil(
        recip_lattice, win["mp_grid"], win["kpoints"]; atol = win["kmesh_tol"]
    )

    ref_bvectors = vec3.([
        [0.0, 0.0, 0.180597],
        [0.0, 0.0, -0.180597],
        [0.188176, 0.0, 0.000001],
        [-0.188176, 0.0, -0.000001],
        [-0.094088, 0.162971, -0.0],
        [0.094088, 0.162971, 0.0],
        [0.094088, -0.162971, 0.0],
        [-0.094088, -0.162971, -0.0],
        [-0.188176, 0.0, 0.180597],
        [0.188176, 0.0, -0.180597],
    ])
    ref_bweights = [
        15.330158,
        15.330158,
        9.413752,
        9.413752,
        9.412853,
        9.412853,
        9.412853,
        9.412853,
        0.000048,
        0.000048,
    ]
    ref_kstencil = Wannier.KspaceStencil(
        recip_lattice,
        vec3(win["mp_grid"]),
        win["kpoints"],
        ref_bvectors,
        ref_bweights,
        nnkp.kpb_k,
        nnkp.kpb_G,
    )
    @test isapprox(kstencil, ref_kstencil; atol = 1.0e-6)
end

@testitem "has_cubic_neighbors" begin
    using Wannier.Datasets
    f = dataset"SnSe2/outputs/SnSe2.nnkp"
    @test Wannier.has_cubic_neighbors(f) == true

    f = dataset"CuBr2/outputs/CuBr2.nnkp"
    @test Wannier.has_cubic_neighbors(f) == false
end

@testitem "compute_bweights" begin
    using Wannier: compute_bweights

    # This set of bvectors are not ordered by norm
    bvectors = vec3.([
        [0.0, 0.0, 0.3267379723643249],
        [0.41527557136667026, 0.0, 0.16336898618216245],
        [0.0, -0.41527557136667026, 0.16336898618216245],
        [0.0, 0.41527557136667026, 0.16336898618216245],
        [-0.41527557136667026, 0.0, 0.16336898618216245],
        [0.0, 0.0, -0.3267379723643249],
        [-0.41527557136667026, 0.0, -0.16336898618216245],
        [0.0, 0.41527557136667026, -0.16336898618216245],
        [0.0, -0.41527557136667026, -0.16336898618216245],
        [0.41527557136667026, 0.0, -0.16336898618216245],
    ])
    weights = compute_bweights(bvectors)

    ref_weights = [
        3.23383985750139,
        1.4496635932248616,
        1.4496635932248616,
        1.4496635932248616,
        1.4496635932248616,
        3.23383985750139,
        1.4496635932248616,
        1.4496635932248616,
        1.4496635932248616,
        1.4496635932248616,
    ]
    @test isapprox(weights, ref_weights; atol = 1.0e-6)
end

@testitem "higher_order_coefficients" begin
    # Σ_m c_m m^(2k) = δ_k1 for k = 1..n, which is what makes the replicated
    # stencil satisfy every even moment up to 2n
    for n in 1:5
        c = Wannier.higher_order_coefficients(n)
        @test length(c) == n
        for k in 1:n
            @test sum(c[m] * m^(2k) for m in 1:n) ≈ (k == 1 ? 1 : 0) atol = 1.0e-12
        end
    end
    # wannier90's values for higher_order_n = 2
    @test Wannier.higher_order_coefficients(2) ≈ [4 / 3, -1 / 12]
end

@testitem "generate_kspace_stencil higher order" begin
    using LinearAlgebra
    using Wannier.Datasets
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = Wannier.reciprocal_lattice(win["unit_cell_cart"])
    gen(; kwargs...) = generate_kspace_stencil(
        recip_lattice, win["mp_grid"], win["kpoints"]; kwargs...
    )

    # `order = 1` must not perturb wannier90's default at all
    ref = gen()
    s1 = gen(; order = 1)
    @test s1.bweights == ref.bweights
    @test s1.kpb_k == ref.kpb_k
    @test s1.kpb_G == ref.kpb_G

    # moments of the stencil: ∑_b w_b b^⊗n
    function moment(stencil, n)
        M = zeros(ntuple(_ -> 3, n))
        for (wb, b) in zip(stencil.bweights, stencil.bvectors)
            for idx in CartesianIndices(M)
                M[idx] += wb * prod(b[i] for i in Tuple(idx))
            end
        end
        return M
    end

    # first order satisfies the B1 condition but not the 4th moment
    @test isapprox(moment(s1, 2), I; atol = 1.0e-6)
    @test maximum(abs, moment(s1, 4)) > 1.0e-3 * norm(s1.bvectors[1])^4

    s2 = gen(; order = 2)
    # the shells are replicated at b and 2b, so the b-vector count doubles
    @test n_bvectors(s2) == 2 * n_bvectors(s1)
    @test isapprox(
        sort(unique(round.(norm.(s2.bvectors); digits = 6))),
        [1, 2] .* norm(s1.bvectors[1]); atol = 1.0e-6,
    )
    # weights are the first-order ones scaled by the analytic coefficients
    @test isapprox(
        sort(unique(round.(s2.bweights; digits = 6))),
        sort(s1.bweights[1] .* Wannier.higher_order_coefficients(2)); atol = 1.0e-6,
    )
    # the 2b shell carries a negative weight, the signature of higher-order FD
    @test any(<(0), s2.bweights)
    # and the 4th moment is annihilated
    @test isapprox(moment(s2, 2), I; atol = 1.0e-6)
    @test maximum(abs, moment(s2, 4)) < 1.0e-6 * norm(s2.bvectors[1])^4

    # third order kills the 6th moment too
    s3 = gen(; order = 3)
    @test n_bvectors(s3) == 3 * n_bvectors(s1)
    @test isapprox(moment(s3, 2), I; atol = 1.0e-6)
    @test maximum(abs, moment(s3, 4)) < 1.0e-6 * norm(s3.bvectors[1])^4
    @test maximum(abs, moment(s3, 6)) < 1.0e-6 * norm(s3.bvectors[1])^6
end

@testitem "compute_bweights higher order from bvectors" begin
    using Wannier.Datasets
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = Wannier.reciprocal_lattice(win["unit_cell_cart"])

    # the `mmn` path only has the b-vector list, so the weights must be
    # recovered from it. With shells at |b| and 2|b| the first-order condition
    # alone is underdetermined, so the higher moments are what pin them down.
    for order in (2, 3)
        s = generate_kspace_stencil(
            recip_lattice, win["mp_grid"], win["kpoints"]; order
        )
        recovered = Wannier.KspaceStencil(
            recip_lattice, win["kpoints"], s.kpb_k, s.kpb_G; order
        )
        @test isapprox(recovered.bweights, s.bweights; rtol = 1.0e-6)
    end

    # recovering with the wrong order is ill-posed and must fail loudly rather
    # than return the minimum-norm solution
    s2 = generate_kspace_stencil(
        recip_lattice, win["mp_grid"], win["kpoints"]; order = 2
    )
    @test_throws AssertionError Wannier.KspaceStencil(
        recip_lattice, win["kpoints"], s2.kpb_k, s2.kpb_G; order = 1
    )
end
