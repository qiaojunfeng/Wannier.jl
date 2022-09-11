
model_ws = read_w90_post(joinpath(FIXTURE_PATH, "valence/band/silicon"); mdrs=false)
model_mdrs = read_w90_post(joinpath(FIXTURE_PATH, "valence/band/silicon"); mdrs=true)

@testset "Hamiltonian WS" begin
    ref_band = WannierIO.read_w90_band(joinpath(FIXTURE_PATH, "valence/band/ws/silicon"))
    E = Wannier.interpolate(model_ws, ref_band.kpoints)
    @test all(isapprox.(E, ref_band.E; atol=2e-5))
end

@testset "Hamiltonian MDRS v2" begin
    ref_band = WannierIO.read_w90_band(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))
    E = Wannier.interpolate(model_mdrs, ref_band.kpoints)
    @test all(isapprox.(E, ref_band.E; atol=2e-5))
end
