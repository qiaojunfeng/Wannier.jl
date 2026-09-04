@testitem "InterpolationModel Wigner-Seitz bands" begin
    using Wannier.Datasets

    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/outputs/Si2_valence.chk.fmt"
    )
    interpolation_model = InterpolationModel(model; real_space = WignerSeitz())

    @test keys(interpolation_model.operators) == (:hamiltonian,)
    @test issorted(interpolation_model.real_space.vectors; by = Tuple)
    @test length(interpolation_model.real_space.vector_index) == n_Rvectors(interpolation_model)
    @test size(interpolation_model.operators.hamiltonian.coefficients) ==
        (n_wannier(model), n_wannier(model), n_Rvectors(interpolation_model))
    @test interpolation_model.operators.hamiltonian.hermitian

    hamiltonian_on_mesh = Wannier._evaluate_real_space_operator(
        interpolation_model.operators.hamiltonian,
        interpolation_model.real_space,
        model.kstencil.kpoints,
        eachindex(model.kstencil.kpoints),
    )
    # The text win file stores the fractional mesh coordinates to eight digits.
    @test maximum(
        abs,
        hamiltonian_on_mesh -
            Wannier.transform_gauge(model.eigenvalues, model.gauges),
    ) < 5.0e-7

    reference = read_w90_band(dataset"Si2_valence/outputs/WS/Si2_valence")
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    segment = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    path = KPath(segment)
    result = interpolate(interpolation_model, path.points, BandEnergy())
    @test keys(result) == (:band_energy,)
    @test all(
        isapprox(view(result.band_energy, :, index), reference.eigenvalues[index]; atol = 2.0e-6)
            for index in axes(result.band_energy, 2)
    )

    bands = band_structure(interpolation_model, path)
    @test bands.band_energy == result.band_energy
    @test bands.kpoints == path.points
    @test length(bands.path_coordinate) == length(path)
end

@testitem "InterpolationModel minimum-distance bands and tensor storage" begin
    using Wannier.Datasets
    using LinearAlgebra: Diagonal

    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/outputs/Si2_valence.chk.fmt"
    )
    number_bands = n_bands(model)
    number_kpoints = n_kpoints(model)
    vector_values = zeros(ComplexF64, number_bands, number_bands, 3, number_kpoints)
    for kpoint_index in 1:number_kpoints, component_index in 1:3
        vector_values[:, :, component_index, kpoint_index] .=
            Diagonal(fill(component_index, number_bands))
    end
    test_vector = BlochOperator(
        vector_values; law = PolarVector(time_reversal = Even()), hermitian = true
    )
    interpolation_model = InterpolationModel(
        model; operators = (; test_vector), real_space = MinimumDistance()
    )

    @test keys(interpolation_model.operators) == (:hamiltonian, :test_vector)
    @test component_shape(interpolation_model.operators.test_vector) == (3,)
    @test interpolation_model.operators.test_vector.hermitian
    @test size(interpolation_model.operators.test_vector.coefficients) == (
        n_wannier(model), n_wannier(model), 3, n_Rvectors(interpolation_model),
    )
    @test all(
        n_Rvectors(operator) == n_Rvectors(interpolation_model)
            for operator in values(interpolation_model.operators)
    )

    evaluated_vector = Wannier._evaluate_real_space_operator(
        interpolation_model.operators.test_vector,
        interpolation_model.real_space,
        model.kstencil.kpoints,
        1:3,
    )
    for kpoint_index in 1:3, component_index in 1:3
        @test isapprox(
            evaluated_vector[:, :, component_index, kpoint_index],
            Diagonal(fill(component_index, n_wannier(model)));
            atol = 2.0e-6,
        )
    end

    reference = read_w90_band(dataset"Si2_valence/outputs/MDRS/Si2_valence")
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    segment = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    path = KPath(segment)
    result = interpolate(interpolation_model, path.points, (BandEnergy(),))
    @test all(
        isapprox(view(result.band_energy, :, index), reference.eigenvalues[index]; atol = 1.0e-7)
            for index in axes(result.band_energy, 2)
    )
end

@testitem "Wannier90 file-data Hamiltonian adapters" begin
    using Wannier.Datasets

    source_model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence",
        dataset"Si2_valence/outputs/Si2_valence.chk.fmt",
    )
    inverse_lattice = inv(source_model.lattice)
    fractional_centers = map(Wannier.center(source_model)) do cartesian_center
        inverse_lattice * cartesian_center
    end
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    segment = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    path = KPath(segment)

    for (directory, tolerance) in (("WS", 2.0e-6), ("MDRS", 2.0e-7))
        prefix = joinpath(
            dataset"Si2_valence/outputs", directory, "Si2_valence"
        )
        tbdat = read_w90_tb_dat(prefix * "_tb.dat")
        hrdat = read_w90_hr_dat(prefix * "_hr.dat")
        wsvec = read_w90_wsvec_dat(prefix * "_wsvec.dat")

        from_tb = InterpolationModel(tbdat; fractional_centers, wsvec)
        from_hr = InterpolationModel(
            hrdat, source_model.lattice; fractional_centers, wsvec
        )
        @test from_tb.real_space.vectors == from_hr.real_space.vectors
        @test maximum(
            abs,
            from_tb.operators.hamiltonian.coefficients -
                from_hr.operators.hamiltonian.coefficients,
        ) < 6.0e-7

        reference = read_w90_band(prefix)
        result_tb = interpolate(from_tb, path.points, BandEnergy())
        result_hr = interpolate(from_hr, path.points, BandEnergy())
        @test all(
            isapprox(
                    view(result_tb.band_energy, :, index),
                    reference.eigenvalues[index];
                    atol = tolerance,
                )
                for index in axes(result_tb.band_energy, 2)
        )
        # The formatted hr.dat stores fewer decimal places than tb.dat.
        @test all(
            isapprox(
                    view(result_hr.band_energy, :, index),
                    reference.eigenvalues[index];
                    atol = 5.0e-5,
                )
                for index in axes(result_hr.band_energy, 2)
        )
    end
end
