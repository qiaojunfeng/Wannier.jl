
@testset "read_w90" begin
    model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))

    @test model.n_bands ≈ 12
    @test model.n_wann ≈ 8
    @test model.n_kpts ≈ 64
end

@testset "write_w90" begin
    model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))

    outdir = mktempdir(; cleanup=true)
    outseedname = joinpath(outdir, "silicon")
    write_w90(outseedname, model)

    A = read_orthonorm_amn("$outseedname.amn")
    @test A ≈ model.A

    E = read_eig("$outseedname.eig")
    @test E ≈ model.E

    M, kpb_k, kpb_b = read_mmn("$outseedname.mmn")
    @test M ≈ model.M
    @test kpb_k ≈ model.bvectors.kpb_k
    @test kpb_b ≈ model.bvectors.kpb_b
end

@testset "read_w90_post" begin
    model = read_w90_post(joinpath(FIXTURE_PATH, "valence/band/silicon"))

    @test typeof(model) <: Wannier.InterpolationModel
    @test model.model.n_wann == 4
    @test model.model.n_kpts == 64
end
