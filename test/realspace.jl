using Printf

const FIXTURE_PATH = joinpath(dirname(pathof(Wannier)), "..", "test", "fixtures")

# A reusable fixture for a model
model = read_w90(joinpath(FIXTURE_PATH, "graphene_unk/graphene"))
model.gauges .= read_amn(joinpath(FIXTURE_PATH, "graphene_unk/graphene.w90.amn"))
unkdir = joinpath(FIXTURE_PATH, "graphene_unk")

@testitem "realspace xsf" begin
    outdir = mktempdir(; cleanup=true)
    outseedname = joinpath(outdir, "wjl")

    write_realspace_wf(outseedname, model; n_supercells=2, unkdir=unkdir)

    for i in 1:(model.n_wann)
        outxsf = read_xsf(joinpath(outdir, @sprintf("wjl_%05d.xsf", i)))
        refxsf = read_xsf(joinpath(unkdir, @sprintf("wjl_%05d.xsf", i)))

        @test isapprox(outxsf.W, refxsf.W; atol=1e-5)
        @test isapprox(outxsf.atom_positions, refxsf.atom_positions; atol=1e-5)
        @test outxsf.atoms == refxsf.atoms
        @test isapprox(outxsf.convvec, refxsf.convvec; atol=1e-5)
        @test isapprox(outxsf.primvec, refxsf.primvec; atol=1e-5)
        @test isapprox(outxsf.rgrid.X, refxsf.rgrid.X; atol=1e-5)
        @test isapprox(outxsf.rgrid.Y, refxsf.rgrid.Y; atol=1e-5)
        @test isapprox(outxsf.rgrid.Z, refxsf.rgrid.Z; atol=1e-5)
        @test isapprox(outxsf.rgrid.basis, refxsf.rgrid.basis; atol=1e-5)
    end
end

@testitem "wannier function" begin
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

@testitem "realspace moment" begin
    rgrid, W = Wannier.read_realspace_wf(model, 2, unkdir)
    r = center.(W)
    ref_r = Vec3[
        Vec3(-0.06333672367474275, -0.0759789055814008, -0.45454804121365266),
        Vec3(0.45683883364060957, 0.49449360532522013, -0.6334692423594555),
        Vec3(0.2941315329655411, 0.49055379163535956, 0.21573045041594663),
        Vec3(0.9742717088741712, 0.006146801998552827, 0.21856691653369997),
        Vec3(0.30104954709081944, 0.46145949408356657, 0.21606232486371238),
    ]
    @test all(isapprox.(r, ref_r; atol=1e-4))

    ω = omega.(W)
    ref_ω = [
        26.87076039248946,
        28.335711581587734,
        25.84537354369948,
        29.79094574412249,
        28.9947646543186,
    ]
    @test isapprox(ω, ref_ω; atol=1e-3)
end
