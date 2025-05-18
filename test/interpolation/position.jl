@testitem "Position WS" begin
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/outputs/Si2_valence.chk.fmt"
    )
    Rspace = generate_Rspace(model; MDRS=false)
    position = TBPosition(Rspace, model; force_hermiticity=false)

    _, ref_position = read_w90_tb(dataset"Si2_valence/outputs/WS/Si2_valence")
    @test all(isapprox.(position, ref_position; atol=2e-8))
end
