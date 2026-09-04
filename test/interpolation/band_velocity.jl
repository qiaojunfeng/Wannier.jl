@testitem "BandEnergy and BandVelocity share one plan" begin
    using Wannier.Datasets

    model = load_dataset("Si2_valence")
    interpolation_model = InterpolationModel(
        model; real_space = MinimumDistance()
    )
    kpoints = [Vec3(0.1, 0.2, 0.3), Vec3(0.17, -0.08, 0.23)]

    energy_only = interpolate(interpolation_model, kpoints, BandEnergy())
    velocity_only = interpolate(interpolation_model, kpoints, BandVelocity())
    combined = interpolate(
        interpolation_model, kpoints, (BandEnergy(), BandVelocity())
    )

    @test keys(combined) == (:band_energy, :band_velocity)
    @test size(combined.band_energy) == (n_wannier(model), length(kpoints))
    @test size(combined.band_velocity) == (3, n_wannier(model), length(kpoints))
    @test combined.band_energy == energy_only.band_energy
    @test combined.band_velocity == velocity_only.band_velocity

    reciprocal = reciprocal_lattice(model)
    step = 1.0e-5
    for cartesian_component in 1:3
        direction = Vec3(inv(reciprocal)[:, cartesian_component])
        forward_kpoints = map(kpoint -> kpoint + step * direction, kpoints)
        backward_kpoints = map(kpoint -> kpoint - step * direction, kpoints)
        forward = interpolate(interpolation_model, forward_kpoints, BandEnergy())
        backward = interpolate(interpolation_model, backward_kpoints, BandEnergy())
        finite_difference = (forward.band_energy - backward.band_energy) / (2step)
        @test isapprox(
            combined.band_velocity[cartesian_component, :, :],
            finite_difference;
            atol = 2.0e-8,
        )
    end

    plan = Wannier._plan_interpolation(
        interpolation_model, (BandEnergy(), BandVelocity()), length(kpoints)
    )
    workspace = Wannier._allocate_interpolation_workspace(plan)
    destination = Wannier._allocate_interpolation_result(
        interpolation_model, plan.observables, length(kpoints)
    )
    Wannier.interpolate!(destination, plan, kpoints, workspace)
    @test destination == combined

    energy_plan = Wannier._plan_interpolation(
        interpolation_model, (BandEnergy(),), length(kpoints)
    )
    energy_workspace = Wannier._allocate_interpolation_workspace(energy_plan)
    @test isempty(keys(energy_workspace.primitive_derivatives))
    @test keys(workspace.primitive_derivatives) == (:hamiltonian,)

    @test_throws ArgumentError interpolate(
        interpolation_model, kpoints, (BandEnergy(), BandEnergy())
    )
    @test_throws ArgumentError interpolate(interpolation_model, kpoints, ())

    no_hamiltonian = InterpolationModel(
        interpolation_model.crystal,
        interpolation_model.basis,
        interpolation_model.real_space,
        (; test_operator = interpolation_model.operators.hamiltonian),
        nothing,
    )
    missing_operator_error = try
        interpolate(no_hamiltonian, kpoints, BandEnergy())
        nothing
    catch error
        error
    end
    @test missing_operator_error isa ArgumentError
    @test occursin("BandEnergy", sprint(showerror, missing_operator_error))
    @test occursin(":hamiltonian", sprint(showerror, missing_operator_error))

    function steady_state_allocations(number_kpoints)
        repeated_kpoints = fill(first(kpoints), number_kpoints)
        allocation_plan = Wannier._plan_interpolation(
            interpolation_model,
            (BandEnergy(), BandVelocity()),
            number_kpoints,
        )
        allocation_workspace = Wannier._allocate_interpolation_workspace(
            allocation_plan
        )
        allocation_result = Wannier._allocate_interpolation_result(
            interpolation_model, allocation_plan.observables, number_kpoints
        )
        Wannier.interpolate!(
            allocation_result,
            allocation_plan,
            repeated_kpoints,
            allocation_workspace,
        )
        return @allocated Wannier.interpolate!(
            allocation_result,
            allocation_plan,
            repeated_kpoints,
            allocation_workspace,
        )
    end

    # Fixed per-batch view bookkeeping is allowed; no allocation scales per k point.
    @test steady_state_allocations(16) <= steady_state_allocations(1) + 256
end
