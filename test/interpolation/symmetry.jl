@testitem "symmetry-closed real-space interpolation" begin
    using WannierIO, Wannier.Datasets
    using LinearAlgebra, Random

    nnkp = read_nnkp(dataset"Si2_hse/outputs/Si2.nnkp")
    stencil = Wannier.KSpaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym(dataset"Si2_hse/Si2.isym")
    Wannier.rescale_littlegroup_reps!(isym.littlegroup_reps)
    eigenvalues_ibz = read_eig(dataset"Si2_hse/Si2.ieig")
    clean_littlegroup_reps!(isym.littlegroup_reps, eigenvalues_ibz)
    centers = [projection.center for projection in nnkp["projections"]]
    constraint = SymmetryConstraint(stencil, isym, centers)

    global_stencil = Wannier.globalize_bvector_ordering(stencil)
    overlaps_ibz = read_mmn(dataset"Si2_hse/Si2.immn").M
    overlaps_fbz = Wannier.reconstruct_overlaps(overlaps_ibz, constraint)
    projections_ibz = read_amn(dataset"Si2_hse/Si2.iamn").A
    gauge_ibz = Wannier.project_covariant(projections_ibz, constraint)
    gauge_fbz = Wannier.reconstruct_gauges(gauge_ibz, constraint)
    eigenvalues_fbz = Wannier.unfold_eigvals(
        eigenvalues_ibz, [collect(pair) for pair in constraint.fbz2ibz]
    )
    win = read_win(dataset"Si2_hse/Si2.win")
    frozen = Wannier.get_frozen_bands(
        eigenvalues_fbz, get(win, "dis_froz_max", -Inf)
    )
    atom_positions = [atom.second for atom in win["atoms_frac"]]
    atom_labels = [string(atom.first) for atom in win["atoms_frac"]]
    model = Wannier.Model(
        win["unit_cell_cart"],
        atom_positions,
        atom_labels,
        global_stencil,
        overlaps_fbz,
        gauge_fbz,
        eigenvalues_fbz,
        frozen,
    )

    raw_model = InterpolationModel(model; real_space = WignerSeitz())
    closed_model = InterpolationModel(
        model;
        real_space = WignerSeitz(),
        symmetry = constraint.wannier_symmetry,
    )
    @test closed_model.symmetry === constraint.wannier_symmetry
    @test closed_model.basis.fractional_centers == centers
    @test Wannier.n_Rvectors(closed_model) > Wannier.n_Rvectors(raw_model)

    kpoints = [Vec3(0.137, -0.083, 0.219), Vec3(-0.071, 0.191, 0.043)]
    for kpoint in kpoints
        reference = Wannier._project_operator_at_kpoint(
            raw_model.operators.hamiltonian,
            raw_model.real_space,
            model.lattice,
            constraint.wannier_symmetry,
            kpoint,
        )
        evaluated = view(
            Wannier._evaluate_real_space_operator(
                closed_model.operators.hamiltonian,
                closed_model.real_space,
                [kpoint],
                Base.OneTo(1),
            ),
            :, :, 1,
        )
        @test norm(evaluated - reference) < 5.0e-12
        @test norm(evaluated - evaluated') < 1.0e-12
    end
    diagnostic = Wannier._symmetry_covariance_residual(
        closed_model, :hamiltonian, kpoints
    )
    @test diagnostic.maximum < 1.0e-11
    @test diagnostic.kpoint_index in eachindex(kpoints)
    @test diagnostic.symmetry_index in eachindex(isym.symops)

    sample_indices = 1:12
    mesh_hamiltonian = Wannier._evaluate_real_space_operator(
        closed_model.operators.hamiltonian,
        closed_model.real_space,
        model.kstencil.kpoints,
        sample_indices,
    )
    mesh_reference = Wannier.transform_gauge(model.eigenvalues, model.gauges)
    # The projection changes the sampled matrices only at the cleaned .isym
    # source-data noise floor.
    @test maximum(
        abs, mesh_hamiltonian - mesh_reference[:, :, sample_indices]
    ) < 2.0e-8

    gamma_hamiltonian = view(
        Wannier._evaluate_real_space_operator(
            closed_model.operators.hamiltonian,
            closed_model.real_space,
            [Vec3(0.0, 0.0, 0.0)],
            Base.OneTo(1),
        ),
        :, :, 1,
    )
    gamma_energies = eigvals(Hermitian(gamma_hamiltonian))
    # The two cubic triplets are algebraically degenerate after closure.
    @test maximum(abs, diff(gamma_energies[2:4])) < 2.0e-14
    @test maximum(abs, diff(gamma_energies[6:8])) < 2.0e-14

    # A small magnetic subgroup exercises nontrivial spatial rotations and
    # antiunitary parity without making every component-law test traverse the
    # full Si2 symmetry group.
    source_indices = [1, 25, 49, 73] # identity, inversion, T, inversion*T
    subgroup_operations = map(enumerate(source_indices)) do (new_index, source_index)
        source = isym.symops[source_index]
        WannierIO.SymOp(
            source.comment,
            source.W,
            source.v,
            source.Wk,
            source.time_reversal,
            source.u,
            new_index,
            new_index,
        )
    end
    number_wannier = length(centers)
    subgroup_representations = map(enumerate(source_indices)) do (new_index, source_index)
        WannierIO.OrbitalRep{number_wannier}(
            new_index, isym.orbital_reps[source_index].D
        )
    end
    subgroup = WannierSymmetry(subgroup_operations, subgroup_representations, centers)
    vectors = [Vec3(-1, 0, 0), Vec3(0, 0, 0), Vec3(1, 0, 0)]
    domain = Wannier.RealSpaceDomain(model.lattice, vectors)
    Random.seed!(1234)

    vector_values = randn(
        ComplexF64,
        Wannier.n_bands(model),
        Wannier.n_bands(model),
        3,
        Wannier.n_kpoints(model),
    )
    vector_operator = BlochOperator(
        vector_values;
        law = PolarVector(time_reversal = Even()),
        hermitian = false,
    )
    multi_operator_model = InterpolationModel(
        model;
        operators = (; test_vector = vector_operator),
        real_space = WignerSeitz(),
        symmetry = subgroup,
    )
    @test keys(multi_operator_model.operators) == (:hamiltonian, :test_vector)
    @test all(
        Wannier.n_Rvectors(operator) == Wannier.n_Rvectors(multi_operator_model)
            for operator in values(multi_operator_model.operators)
    )
    @test Wannier._symmetry_covariance_residual(
        multi_operator_model, :test_vector, kpoints
    ).maximum < 1.0e-11

    for law in (
            PolarVector(time_reversal = Even()),
            AxialVector(time_reversal = Odd()),
            CartesianTensor(2; time_reversal = Even()),
        )
        shape = component_shape(law)
        coefficients = randn(
            ComplexF64, number_wannier, number_wannier, shape..., length(vectors)
        )
        raw_operator = Wannier.RealSpaceOperator(
            coefficients, law, domain; hermitian = false
        )
        realization = Wannier._close_real_space_operator(
            vectors,
            (; law, hermitian = false),
            coefficients,
            subgroup,
            model.lattice,
        )
        closed_vectors = sort!(collect(keys(realization)); by = Tuple)
        component_domain = Wannier.RealSpaceDomain(model.lattice, closed_vectors)
        closed_coefficients = Wannier._pack_real_space_dictionary(
            realization, closed_vectors, law
        )
        closed_operator = Wannier.RealSpaceOperator(
            closed_coefficients, law, component_domain; hermitian = false
        )
        component_model = InterpolationModel(
            closed_model.crystal,
            closed_model.basis,
            component_domain,
            (; test_operator = closed_operator),
            subgroup,
        )

        reference = Wannier._project_operator_at_kpoint(
            raw_operator, domain, model.lattice, subgroup, first(kpoints)
        )
        evaluated = Wannier._evaluate_real_space_operator(
            closed_operator, component_domain, [first(kpoints)], Base.OneTo(1)
        )
        evaluated = dropdims(evaluated; dims = ndims(evaluated))
        @test norm(evaluated - reference) < 5.0e-12
        @test Wannier._symmetry_covariance_residual(
            component_model, :test_operator, kpoints
        ).maximum < 1.0e-11
    end

    # Exercise tied minimum-distance replicas before symmetry closure. All
    # copies of one (row, column, quotient-vector) contribution retain the
    # same scheme-assigned weight, and the closed Fourier sum still equals the
    # independent query projector.
    selection = Wannier._real_space_selection(
        MinimumDistance(), model.lattice, (1, 1, 1), centers
    )
    quotient_matrix = randn(ComplexF64, number_wannier, number_wannier)
    quotient_matrix = (quotient_matrix + quotient_matrix') / 2
    quotient_coefficients = reshape(quotient_matrix, number_wannier, number_wannier, 1)
    selected_coefficients = Wannier._apply_real_space_selection(
        selection, quotient_coefficients
    )
    replica_groups = Dict{Tuple{Int, Int, Int}, Vector{Tuple{Int, Int}}}()
    for (representative_index, contributions) in enumerate(selection.contributions)
        for contribution in contributions
            key = (
                contribution.quotient_index, contribution.row, contribution.column,
            )
            push!(
                get!(replica_groups, key, Tuple{Int, Int}[]), (
                    representative_index, contribution.denominator,
                )
            )
        end
    end
    matching_key = first(
        key for (key, group) in replica_groups if length(group) > 1
    )
    tied_group = replica_groups[matching_key]
    @test length(tied_group) > 1
    _, denominator = first(tied_group)
    _, row, column = matching_key
    expected_weighted_value = quotient_matrix[row, column] / denominator
    @test all(tied_group) do (representative_index, replica_denominator)
        replica_denominator == denominator &&
            selected_coefficients[row, column, representative_index] ≈
            expected_weighted_value
    end

    selected_vectors = Wannier._representative_vectors(selection)
    selected_domain = Wannier.RealSpaceDomain(model.lattice, selected_vectors)
    selected_order = Wannier._canonical_vector_order(selected_domain, selected_vectors)
    selected_coefficients_ordered = Wannier._reorder_vector_axis(
        selected_coefficients, selected_order
    )
    selected_operator = Wannier.RealSpaceOperator(
        selected_coefficients_ordered,
        Scalar(time_reversal = Even()),
        selected_domain;
        hermitian = false,
    )
    tied_realization = Wannier._close_real_space_operator(
        selected_vectors,
        (; law = selected_operator.law, hermitian = false),
        selected_coefficients,
        subgroup,
        model.lattice,
    )
    tied_vectors = sort!(collect(keys(tied_realization)); by = Tuple)
    tied_domain = Wannier.RealSpaceDomain(model.lattice, tied_vectors)
    tied_coefficients = Wannier._pack_real_space_dictionary(
        tied_realization, tied_vectors, selected_operator.law
    )
    tied_operator = Wannier.RealSpaceOperator(
        tied_coefficients, selected_operator.law, tied_domain; hermitian = false
    )
    tied_reference = Wannier._project_operator_at_kpoint(
        selected_operator, selected_domain, model.lattice, subgroup, first(kpoints)
    )
    tied_evaluated = view(
        Wannier._evaluate_real_space_operator(
            tied_operator, tied_domain, [first(kpoints)], Base.OneTo(1)
        ),
        :, :, 1,
    )
    @test norm(tied_evaluated - tied_reference) < 5.0e-12
end
