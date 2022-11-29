
@testset "Fermi surface" begin
    model = read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))

    ref_bxsf = WannierIO.read_bxsf(joinpath(FIXTURE_PATH, "valence/band/mdrs/wjl.bxsf"))

    kpoints, E = Wannier.fermi_surface(model; n_k=2)

    @test all(isapprox.(E, ref_bxsf.E; atol=1e-7))
end
