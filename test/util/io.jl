
@testset "isbinary file" begin
    @test !Wannier.isbinary(joinpath(FIXTURE_PATH, "silicon/UNK00001.1"))
end

@testset "parse_float" begin
    @test Wannier.parse_float("1.0D-10") ≈ 1e-10
end
