
@testset "isbinary_file" begin
    @test !Wannier.isbinary_file("$FIXTURE_PATH/UNK00001.1")
end
