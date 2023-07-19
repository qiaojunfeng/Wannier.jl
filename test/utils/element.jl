@testitem "get_atom_number" begin
    @test Wannier.get_atom_number("O") == 8
    @test Wannier.get_atom_number(["H", "O"]) == [1, 8]
end
