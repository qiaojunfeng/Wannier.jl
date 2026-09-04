@testitem "spin projection" begin
    using LinearAlgebra
    using Wannier.Datasets
    hamiltonian, position, spin = read_w90_tb_chk_spn(
        dataset"Fe_soc_coarse/outputs/Fe";
        spn = dataset"Fe_soc_coarse/Fe.spn",
        chk = dataset"Fe_soc_coarse/outputs/Fe.chk",
    )
    # project onto the z axis
    θ = 0.0
    ϕ = 0.0
    interp = SpinProjectionInterpolator(hamiltonian, spin, θ, ϕ)

    ref_kpt = read_w90_band_kpt(dataset"Fe_soc_coarse/outputs/postw90/Fe-path.kpt")
    ref_dat = read_w90_band_dat(dataset"Fe_soc_coarse/outputs/postw90/Fe-bands.dat")

    # if I use the kpoints in ref_kpt, the difference between eigenvalues is
    # around 1e-4, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_kpt.kpoints
    win = read_win(dataset"Fe_soc_coarse/Fe.win")
    kseg = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    kpath = KPath(kseg, 5)
    # postw90.x has a bug, it misses the `H` point at 22
    kpoints = collect(kpath)
    deleteat!(kpoints, 22)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1.0e-6)
    ##
    eigenvalues = HamiltonianInterpolator(hamiltonian)(kpoints)[1]
    @test all(norm(view(eigenvalues, :, ik) - ref_dat.eigenvalues[ik]) < 2.0e-6 for ik in axes(eigenvalues, 2))

    Sz = interp(kpoints)
    @test all(isapprox.(Sz, ref_dat.extras; atol = 5.0e-5))

    model = read_w90_with_chk(
        dataset"Fe_soc_coarse/Fe", dataset"Fe_soc_coarse/outputs/Fe.chk"
    )
    spin_operator = BlochOperator(read_spn(dataset"Fe_soc_coarse/Fe.spn"))
    interpolation_model = InterpolationModel(
        model;
        operators = (; spin = spin_operator),
        real_space = MinimumDistance(),
    )
    combined = interpolate(
        interpolation_model,
        kpoints,
        (BandEnergy(), SpinExpectation([0.0, 0.0, 1.0]; truncate = true)),
    )
    plan = Wannier._plan_interpolation(
        interpolation_model,
        (BandEnergy(), SpinExpectation([0.0, 0.0, 1.0])),
        length(kpoints),
    )
    workspace = Wannier._allocate_interpolation_workspace(plan)
    @test plan.primitive_operators == (:hamiltonian, :spin)
    @test keys(workspace.primitive_values) == (:spin,)
    @test all(
        norm(view(combined.band_energy, :, index) - ref_dat.eigenvalues[index]) < 2.0e-6
            for index in axes(combined.band_energy, 2)
    )
    @test all(isapprox.(eachcol(combined.spin_expectation), ref_dat.extras; atol = 5.0e-5))

    components = interpolate(interpolation_model, kpoints, SpinExpectation())
    projected = interpolate(
        interpolation_model, kpoints, SpinExpectation([0.0, 0.0, 1.0])
    )
    @test view(components.spin_expectation, 3, :, :) ≈ projected.spin_expectation

    separate = interpolate(
        interpolation_model, kpoints, SpinExpectation([0.0, 0.0, 1.0]; truncate = true)
    )
    @test separate.spin_expectation == combined.spin_expectation
    @test_throws ArgumentError SpinExpectation([0.0, 0.0])
    @test_throws ArgumentError SpinExpectation(zeros(3))
end
