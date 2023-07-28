@testitem "generate_kspace_stencil" begin
    using Wannier: Vec3
    using Wannier.Datasets
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    _, kpb_k, kpb_G = read_mmn(dataset"Si2_valence/Si2_valence.mmn")

    recip_lattice = reciprocal_lattice(win.unit_cell_cart)
    kstencil = generate_kspace_stencil(recip_lattice, win.mp_grid, win.kpoints)

    # copied from wout
    ref_bvectors = Vec3{Float64}[
        [0.192835, 0.192835, -0.192835],
        [0.192835, -0.192835, 0.192835],
        [-0.192835, 0.192835, 0.192835],
        [0.192835, 0.192835, 0.192835],
        [-0.192835, -0.192835, 0.192835],
        [-0.192835, 0.192835, -0.192835],
        [0.192835, -0.192835, -0.192835],
        [-0.192835, -0.192835, -0.192835],
    ]
    ref_bweights = fill(3.361532, 8)
    ref_kstencil = Wannier.KspaceStencil(
        recip_lattice, win.mp_grid, win.kpoints, ref_bvectors, ref_bweights, kpb_k, kpb_G
    )
    @test isapprox(kstencil, ref_kstencil; atol=1e-6)
end

@testitem "generate_kspace_stencil 2D" begin
    using Wannier: Vec3
    using Wannier.Datasets
    win = read_win(dataset"graphene/graphene.win")
    # nnkp = read_nnkp_compute_bweights(dataset"graphene/reference/graphene.nnkp")
    nnkp = read_nnkp_compute_bweights(
        expanduser("~/git/WannierDatasets/datasets/graphene/reference/graphene.nnkp")
    )

    recip_lattice = reciprocal_lattice(win.unit_cell_cart)
    kstencil = generate_kspace_stencil(recip_lattice, win.mp_grid, win.kpoints)

    ref_bvectors = Vec3{Float64}[
        [0.0, 0.0, 0.628319],
        [0.0, 0.0, -0.628319],
        [0.793031, 0.457857, 0.0],
        [-0.793031, -0.457857, 0.0],
        [0.793031, -0.457857, 0.0],
        [-0.793031, 0.457857, 0.0],
        [0.0, -0.915713, 0.0],
        [0.0, 0.915713, 0.0],
    ]
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
        win.mp_grid,
        win.kpoints,
        ref_bvectors,
        ref_bweights,
        nnkp.kpb_k,
        nnkp.kpb_G,
    )
    # wout does not have enough digits, use a bit larger atol
    @test isapprox(kstencil, ref_kstencil; atol=1e-5)
end

@testitem "generate_kspace_stencil kmesh_tol" begin
    using Wannier: Vec3
    using Wannier.Datasets
    win = read_win(expanduser("~/git/WannierDatasets/datasets/SnSe2/SnSe2.win"))
    nnkp = read_nnkp_compute_bweights(
        expanduser("~/git/WannierDatasets/datasets/SnSe2/reference/SnSe2.nnkp")
    )
    recip_lattice = reciprocal_lattice(win.unit_cell_cart)
    kstencil = generate_kspace_stencil(
        recip_lattice, win.mp_grid, win.kpoints; atol=win.kmesh_tol
    )

    ref_bvectors = Vec3{Float64}[
        [0.000000, 0.000000, 0.180597],
        [0.000000, 0.000000, -0.180597],
        [0.188176, 0.000000, 0.000001],
        [-0.188176, 0.000000, -0.000001],
        [-0.094088, 0.162971, -0.000000],
        [0.094088, 0.162971, 0.000000],
        [0.094088, -0.162971, 0.000000],
        [-0.094088, -0.162971, -0.000000],
        [-0.188176, 0.000000, 0.180597],
        [0.188176, 0.000000, -0.180597],
    ]
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
        win.mp_grid,
        win.kpoints,
        ref_bvectors,
        ref_bweights,
        nnkp.kpb_k,
        nnkp.kpb_G,
    )
    @test isapprox(kstencil, ref_kstencil; atol=1e-6)
end
