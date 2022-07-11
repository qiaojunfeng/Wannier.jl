
@testset "cli help" begin
    @test Wannier.command_main(["-h"]) == 0

    for cmd in ["maxloc", "dis", "ptg", "optrot", "splitvc"]
        @test Wannier.command_main([cmd, "-h"]) == 0
    end
end
