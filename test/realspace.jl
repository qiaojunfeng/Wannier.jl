@testmodule RealspaceEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier.Datasets

    model = read_w90_with_chk(dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk")
    unkdir = dataset"Si2_coarse"
end

@testitem "realspace xsf" setup = [RealspaceEnv] begin
    using Printf
    outdir = mktempdir(; cleanup=true)
    outseedname = joinpath(outdir, "wjl")

    write_realspace_wf(outseedname, RealspaceEnv.model; n_supercells=3, RealspaceEnv.unkdir)

    for i in 1:n_wannier(RealspaceEnv.model)
        outxsf = read_xsf(joinpath(outdir, @sprintf("wjl_%05d.xsf", i)))
        refxsf = read_xsf(
            joinpath(RealspaceEnv.unkdir, "outputs", @sprintf("Si2_%05d.xsf", i))
        )

        @test all(isapprox.(outxsf.W, refxsf.W; atol=1e-4))
        @test isapprox(outxsf.atom_positions, refxsf.atom_positions; atol=1e-5)
        # refxsf.atoms = ["si", "si"], written by wannier90
        @test parse.(Int, outxsf.atoms) == Wannier.get_atom_number(titlecase.(refxsf.atoms))
        @test isapprox(outxsf.convvec, refxsf.convvec; atol=1e-5)
        @test isapprox(outxsf.primvec, refxsf.primvec; atol=1e-5)
        @test isapprox(outxsf.rgrid.X, refxsf.rgrid.X; atol=1e-5)
        @test isapprox(outxsf.rgrid.Y, refxsf.rgrid.Y; atol=1e-5)
        @test isapprox(outxsf.rgrid.Z, refxsf.rgrid.Z; atol=1e-5)
        @test isapprox(outxsf.rgrid.basis, refxsf.rgrid.basis; atol=1e-5)
    end
end

@testitem "wannier function" begin
    using LinearAlgebra
    using Wannier: WannierFunction, Vec3, SVector
    x_range = -1:0.01:1
    y_range = -1:0.01:1
    z_range = -1:0.01:1

    wfc_grid = [Vec3(x, y, z) for x in x_range, y in y_range, z in z_range]

    px_orb = normalize(
        Wannier.WannierFunction(
            wfc_grid, [SVector(((p[1] + 0im) * ℯ^(-norm(p)^2)),) for p in wfc_grid]
        ),
    )

    px_orb2 = px_orb'
    LinearAlgebra.adjoint!(px_orb2, px_orb2)

    @test values(px_orb2) == values(px_orb)

    @test norm(px_orb) ≈ norm(px_orb2) ≈ 1.0

    py_orb = normalize(
        WannierFunction(
            wfc_grid, [SVector(((p[2] + 0im) * ℯ^(-norm(p)^2)),) for p in wfc_grid]
        ),
    )

    pz_orb = normalize(
        WannierFunction(
            wfc_grid, [SVector(((p[3] + 0im) * ℯ^(-norm(p)^2)),) for p in wfc_grid]
        ),
    )

    @test dot(px_orb, py_orb) <= 1.0e-15
    @test dot(px_orb, px_orb) ≈ 1.0

    Lx = zeros(ComplexF64, 3, 3)
    Ly = zeros(ComplexF64, 3, 3)
    Lz = zeros(ComplexF64, 3, 3)
    for (i1, p1) in enumerate((px_orb, py_orb, pz_orb)),
        (i2, p2) in enumerate((px_orb, py_orb, pz_orb))

        Lx[i1, i2], Ly[i1, i2], Lz[i1, i2] = Wannier.calc_angmom(p1, p2, zero(Vec3))
    end

    @test norm(sum(Lx .- [0 0 0; 0 0 -im; 0 im 0])) < 1e-4
    @test norm(sum(Ly .- [0 0 im; 0 0 0; -im 0 0])) < 1e-4
    @test norm(sum(Lz .- [0 -im 0; im 0 0; 0 0 0])) < 1e-4

    px_orb_up = normalize(
        WannierFunction(
            wfc_grid,
            [SVector((p[1] + 0im, zero(ComplexF64)) .* ℯ^(-norm(p)^2)) for p in wfc_grid],
        ),
    )
    px_orb_dn = normalize(
        WannierFunction(
            wfc_grid,
            [SVector((zero(ComplexF64), p[1] + 0im) .* ℯ^(-norm(p)^2)) for p in wfc_grid],
        ),
    )

    @test dot(px_orb_dn, px_orb_up) ≈ 0.0
    @test Wannier.calc_spin(px_orb_up, px_orb_up) ≈
        Vec3(0.0 + 0im, 0.0 + 0.0im, 0.5 + 0.0im)

    @test norm(Wannier.calc_dipole(px_orb, py_orb)) < 1e-17
    @test norm(Wannier.calc_dipole(px_orb_up, px_orb_dn)) < 1e-17
end

@testitem "realspace moment" setup = [RealspaceEnv] begin
    using Wannier: Vec3
    rgrid, W = Wannier.read_realspace_wf(RealspaceEnv.model, 2, RealspaceEnv.unkdir)
    Wannier.normalize!.(W)
    r = center.(W)
    ref_r = Vec3[
        [0.05649182474139805, 0.05647898672201181, 0.05648648084336603],
        [0.006835975305507893, 0.006835714628599313, 0.00683543106836339],
        [0.0068280451636423736, 0.0068356659075035275, 0.00683581178070572],
        [0.006835776140180217, 0.006845785639297941, 0.006835766959433606],
        [1.2086580389311699, 1.2086465413956256, 1.2086532094932414],
        [1.0968514169399275, 1.0968512556524777, 1.0968511402797747],
        [1.096844404431083, 1.0968513029064058, 1.0968513855922837],
        [1.0968516077683101, 1.09686041159614, 1.0968514341399862],
    ]
    @test all(isapprox.(r, ref_r; atol=1e-4))

    ω = omega.(W)
    ref_ω = [
        1.925692300941864,
        3.0882234042095624,
        3.0882244433738033,
        3.08822217509327,
        2.488361186325519,
        4.170739171655333,
        4.170733540871117,
        4.170747304039493,
    ]
    @test isapprox(ω, ref_ω; atol=1e-3)
end
