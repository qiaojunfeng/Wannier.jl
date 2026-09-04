module BenchInterpolation

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2_valence")
    kpoints = Wannier.kpoints(model)

    wigner_seitz_model = InterpolationModel(model; real_space = WignerSeitz())
    minimum_distance_model = InterpolationModel(model; real_space = MinimumDistance())

    # Warm up the complete public path before BenchmarkTools measures it.
    interpolate(wigner_seitz_model, kpoints, BandEnergy())
    interpolate(minimum_distance_model, kpoints, BandEnergy())
    interpolate(
        minimum_distance_model, kpoints, (BandEnergy(), BandVelocity())
    )

    path_kpoints = [
        Vec3(t, 0.5t, 0.25t) for t in range(-0.5, 0.5; length = 256)
    ]
    dense_kpoints = [
        Vec3(i / 32, j / 32, k / 32) for k in 0:31 for j in 0:31 for i in 0:31
    ]
    interpolate(minimum_distance_model, path_kpoints, BandEnergy())
    interpolate(minimum_distance_model, dense_kpoints, BandEnergy())

    SUITE["band energy"]["Wigner-Seitz"] =
        @benchmarkable interpolate($wigner_seitz_model, $kpoints, BandEnergy())
    SUITE["band energy"]["minimum distance"] =
        @benchmarkable interpolate($minimum_distance_model, $kpoints, BandEnergy())
    SUITE["band energy"]["path 256 minimum distance"] =
        @benchmarkable interpolate($minimum_distance_model, $path_kpoints, BandEnergy())
    SUITE["band energy"]["dense 32^3 minimum distance"] =
        @benchmarkable interpolate($minimum_distance_model, $dense_kpoints, BandEnergy())

    SUITE["band energy and velocity"]["combined"] = @benchmarkable interpolate(
        $minimum_distance_model,
        $kpoints,
        (BandEnergy(), BandVelocity()),
    )
    SUITE["band energy and velocity"]["separate"] = @benchmarkable begin
        interpolate($minimum_distance_model, $kpoints, BandEnergy())
        interpolate($minimum_distance_model, $kpoints, BandVelocity())
    end

    stage_plan = Wannier._plan_interpolation(
        minimum_distance_model, (BandEnergy(),), length(kpoints)
    )
    stage_workspace = Wannier._allocate_interpolation_workspace(stage_plan)
    stage_indices = 1:length(kpoints)
    stage_phase = view(stage_workspace.phase, :, stage_indices)
    stage_hamiltonian = view(stage_workspace.hamiltonian, :, :, stage_indices)
    stage_eigenvalues = view(stage_workspace.eigenvalues, :, stage_indices)
    stage_eigenvectors = view(stage_workspace.eigenvectors, :, :, stage_indices)
    Wannier._fourier_phase_block!(
        stage_phase, minimum_distance_model.real_space, kpoints
    )
    Wannier._evaluate_real_space_operator!(
        stage_hamiltonian,
        minimum_distance_model.operators.hamiltonian,
        stage_phase,
    )

    SUITE["band-energy stages"]["phase"] = @benchmarkable Wannier._fourier_phase_block!(
        $stage_phase, $(minimum_distance_model.real_space), $kpoints
    )
    SUITE["band-energy stages"]["Hamiltonian Fourier sum"] =
        @benchmarkable Wannier._evaluate_real_space_operator!(
        $stage_hamiltonian,
        $(minimum_distance_model.operators.hamiltonian),
        $stage_phase,
    )
    SUITE["band-energy stages"]["diagonalization"] =
        @benchmarkable Wannier._diagonalize_hermitian_batch!(
        $stage_eigenvalues,
        $stage_eigenvectors,
        $stage_hamiltonian,
        $(stage_workspace.diagonalization),
    )

    fe_prefix = dataset"Fe_soc_coarse/outputs/Fe"
    fe_model = read_w90_tb_chk_spn(
        fe_prefix;
        chk = dataset"Fe_soc_coarse/outputs/Fe.chk",
        spn = dataset"Fe_soc_coarse/Fe.spn",
    )
    fe_fermi_energy = read_win(dataset"Fe_soc_coarse/Fe.win")["fermi_energy"]
    fe_kpoints = [Vec3(t, 0.2t, 0.1) for t in range(-0.5, 0.5; length = 128)]
    fe_observables = (
        BandEnergy(),
        BerryCurvature(fe_fermi_energy),
        SpinExpectation(),
    )
    interpolate(fe_model, fe_kpoints, fe_observables)

    SUITE["energy, Berry curvature, and spin"]["combined"] =
        @benchmarkable interpolate($fe_model, $fe_kpoints, $fe_observables)
    SUITE["energy, Berry curvature, and spin"]["separate"] = @benchmarkable begin
        interpolate($fe_model, $fe_kpoints, BandEnergy())
        interpolate($fe_model, $fe_kpoints, BerryCurvature($fe_fermi_energy))
        interpolate($fe_model, $fe_kpoints, SpinExpectation())
    end

end  # module

BenchInterpolation.SUITE
