using Printf

const FIXTURE_PATH = joinpath(dirname(pathof(Wannier)), "..", "test", "fixtures")

# A reusable fixture for a model
model = read_w90(joinpath(FIXTURE_PATH, "graphene_unk/graphene"))
model.U .= read_amn(joinpath(FIXTURE_PATH, "graphene_unk/graphene.w90.amn"))
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

@testitem "realspace moment" begin
    rgrid, W = read_realspace_wf(model, 2, unkdir)

    r = center(rgrid, W)
    ref_r = [
        -0.13452 0.871192 0.559132 1.96598 0.612977
        -0.16137 0.942999 0.932522 0.0124036 0.939593
        -0.965406 -1.20803 0.410094 0.441045 0.439932
    ]
    ref_r = map(i -> Vec3(ref_r[:, i]), axes(ref_r, 2))
    @test isapprox(r, ref_r; atol=1e-4)

    ω = omega(rgrid, W)
    ref_ω = [56.5538 52.5582 48.4909 58.0671 58.2983]'
    @test isapprox(ω, ref_ω; atol=1e-3)
end
