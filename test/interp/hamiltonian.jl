H_ws = Wannier.read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/ws/silicon")).H;
H_mdrs = Wannier.read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon")).H

# TODO this fails
# @testset "Hamiltonian WS" begin
#     ref_band = WannierIO.read_w90_band(joinpath(FIXTURE_PATH, "valence/band/ws/silicon"))
#     eigenvalues = Wannier.interpolate(H_ws, ref_band.kpoints)
#     @test all(isapprox.(eigenvalues, ref_band.eigenvalues; atol=2e-5))
# end

@testset "Hamiltonian MDRS v2" begin
    ref_band = WannierIO.read_w90_band(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))
    eigenvalues = Wannier.interpolate(H_mdrs, ref_band.kpoints)
    @test all(isapprox.(eigenvalues, ref_band.eigenvalues; atol=2e-5))
end
