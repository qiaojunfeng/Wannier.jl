@testitem "Position WS" begin
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/outputs/Si2_valence.chk.fmt"
    )
    Rspace = generate_Rspace(model; MDRS = false)
    position = TBPosition(Rspace, model; force_hermiticity = false)

    _, ref_position = read_w90_tb(dataset"Si2_valence/outputs/WS/Si2_valence")
    @test all(isapprox.(position, ref_position; atol = 2.0e-8))

    interpolation_model = InterpolationModel(
        model;
        operators = (;
            berry_connection = BerryConnection(; force_hermiticity = false),
        ),
        real_space = WignerSeitz(),
    )
    connection = Wannier._evaluate_real_space_operator(
        interpolation_model.operators.berry_connection,
        interpolation_model.real_space,
        model.kstencil.kpoints,
        eachindex(model.kstencil.kpoints),
    )
    reference_connection = Wannier._dense_berry_connection(
        Wannier.compute_berry_connection_kspace(
            model; force_hermiticity = false
        )
    )
    @test maximum(abs, connection - reference_connection) < 5.0e-7
    @test component_shape(interpolation_model.operators.berry_connection) == (3,)
end
