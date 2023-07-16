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

    U = read_orthonorm_amn("$outseedname.amn")
    @test U ≈ model.U

    E = read_eig("$outseedname.eig")
    @test E ≈ model.E

    M, kpb_k, kpb_G = read_mmn("$outseedname.mmn")
    @test M ≈ model.M
    @test kpb_k ≈ model.bvectors.kpb_k
    @test kpb_G ≈ model.bvectors.kpb_G
end

@testset "read_w90_with_chk" begin
    model = read_w90_with_chk(joinpath(FIXTURE_PATH, "valence/band/silicon"))

    @test model[1] isa Wannier.TBHamiltonian
    @test Wannier.n_wann(model[1]) == 4
end
