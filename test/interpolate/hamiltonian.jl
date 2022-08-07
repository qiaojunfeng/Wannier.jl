
# @testset "interpolate" begin
#     model = read_w90(joinpath(FIXTURE_PATH, "valence/band/silicon"))
#     band = read_w90_band(joinpath(FIXTURE_PATH, "valence/band/silicon"))

#     E = Wannier.interpolate(model, band.kpoints)

#     @test isapprox(E, band.E; atol = 1e-8)
# end
