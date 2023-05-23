model = read_w90(joinpath(FIXTURE_PATH, "valence/band/silicon"); amn=false)
model.U .= get_U(read_chk(joinpath(FIXTURE_PATH, "valence/band/silicon.chk.fmt")))
model_ws = read_w90_tb(
    joinpath(FIXTURE_PATH, "valence/band/ws/silicon"))
model_mdrs = read_w90_tb(
    joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon")
)

@testset "fourier WS" begin
    Hᴿ = Wannier.TBHamiltonian(model, model_ws.R)
    ref_Hᴿ = model_ws.H
    @test all(x-> all(isapprox.(x[1].tb_block, x[2].tb_block;atol=1e-7)), zip(Hᴿ, ref_Hᴿ))
end

# @testset "invfourier WS" begin
#     ref_Hᵏ = Wannier.get_Hk(model.E, model.U)
#     Hᴿ = model_ws.H
#     Hᵏ = map(k -> Wannier.Hk(Hᴿ, k), model.kpoints)
#     Hᵏ[1]
#     ref_Hᵏ[1]
#     @test all(i -> all(isapprox.(Hᵏ[i], ref_Hᵏ[i]; atol=1e-7)), 1:length(ref_Hᵏ))
# end

# @testset "fourier MDRS v1" begin
#     Hᵏ = Wannier.get_Hk(model.E, model.U)
#     Hᴿ = Wannier.fourier(model_mdrs.kRvectors, Hᵏ; version=:v1)
#     ref_Hᴿ = model_ws.H
#     @test all(isapprox.(Hᴿ, ref_Hᴿ; atol=1e-7))
# end

# @testset "invfourier MDRS v1" begin
#     ref_Hᵏ = Wannier.get_Hk(model.E, model.U)
#     Hᴿ = model_ws.H
#     Hᵏ = Wannier.invfourier(model_mdrs.kRvectors, Hᴿ, model.kpoints; version=:v1)
#     @test all(i -> all(isapprox.((@view Hᵏ[:,:,i]), ref_Hᵏ[i]; atol=1e-7)), 1:length(ref_Hᵏ))
# end

# @testset "fourier/invfourier MDRS v2" begin
#     # Hᴿ of MDRS v2 has different R vectors, so I cannot compare with tb.dat
#     ref_Hᵏ = Wannier.get_Hk(model.E, model.U)
#     Hᴿ = Wannier.fourier(model_mdrs.kRvectors, ref_Hᵏ; version=:v2)
#     Hᵏ = Wannier.invfourier(model_mdrs.kRvectors, Hᴿ, model.kpoints; version=:v2)
#     @test all(i -> all(isapprox.((@view Hᵏ[:,:,i]), ref_Hᵏ[i]; atol=1e-7)), 1:length(ref_Hᵏ))
# end
