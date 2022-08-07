model_ws = read_w90_post(joinpath(FIXTURE_PATH, "valence/band/silicon"); mdrs=false)
model_mdrs = read_w90_post(joinpath(FIXTURE_PATH, "valence/band/silicon"); mdrs=true)
tb_ws = read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/ws/silicon"))
tb_mdrs = read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))

# @testset "fourier WS" begin
#     Hᵏ = Wannier.get_Hk(model_ws.model.E, model_ws.model.A)
#     Hᴿ = Wannier.fourier(model_ws.kRvectors, Hᵏ)
#     Rvecs, ref_Hᴿ, positions = tb_ws
#     @test all(isapprox.(Hᴿ, ref_Hᴿ; atol = 1e-7))
# end

# @testset "invfourier WS" begin
#     ref_Hᵏ = Wannier.get_Hk(model_ws.model.E, model_ws.model.A)
#     Rvecs, Hᴿ, positions = tb_ws
#     Hᵏ = Wannier.invfourier(model_ws.kRvectors, Hᴿ, model_ws.model.kpoints)
#     @test all(isapprox.(Hᵏ, ref_Hᵏ; atol = 1e-7))
# end

@testset "fourier MDRS v1" begin
    Hᵏ = Wannier.get_Hk(model_mdrs.model.E, model_mdrs.model.A)
    Hᴿ = Wannier.fourier(model_mdrs.kRvectors, Hᵏ; version=:v1)
    Rvecs, ref_Hᴿ, positions = tb_mdrs
    @test all(isapprox.(Hᴿ, ref_Hᴿ; atol=1e-7))
end

@testset "invfourier MDRS v1" begin
    ref_Hᵏ = Wannier.get_Hk(model_mdrs.model.E, model_mdrs.model.A)
    Rvecs, Hᴿ, positions = tb_mdrs
    Hᵏ = Wannier.invfourier(model_mdrs.kRvectors, Hᴿ, model_mdrs.model.kpoints; version=:v1)
    @test all(isapprox.(Hᵏ, ref_Hᵏ; atol=1e-7))
end
