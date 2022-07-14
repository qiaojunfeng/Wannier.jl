
@testset "isbinary_file" begin
    @test !Wannier.isbinary_file(joinpath(FIXTURE_PATH, "silicon/UNK00001.1"))
end
