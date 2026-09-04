@testitem "BerryCurvature observable" begin
    using LinearAlgebra
    using DelimitedFiles
    using Wannier.Datasets
    model = read_w90_with_chk(dataset"Fe_soc_coarse/Fe", dataset"Fe_soc_coarse/outputs/Fe.chk")
    win = read_win(dataset"Fe_soc_coarse/Fe.win")

    ref_kpt = read_w90_band_kpt(dataset"Fe_soc_coarse/outputs/postw90/Fe-path.kpt")
    ref_dat = readdlm(dataset"Fe_soc_coarse/outputs/postw90/Fe-curv.dat")
    # w90 actually writes -Ω, so we need to negate it
    ref_Ω = map(eachrow(ref_dat[:, 2:end])) do Ω
        -Wannier.axialvector_to_antisymmetrictensor(Ω)
    end

    # if I use the kpoints in ref_kpt, the difference between eigenvalues is
    # around 1e-4, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_kpt.kpoints
    kseg = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    # I use less kpoints to speedup the test
    kpath = KPath(kseg, 5)
    kpoints = collect(kpath)
    # postw90.x has a bug, it misses the `H` point
    ik_H = 22
    deleteat!(kpoints, ik_H)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1.0e-6)

    interpolation_model = InterpolationModel(
        model;
        operators = (;
            berry_connection = BerryConnection(; imlog_diag = false),
        ),
        real_space = MinimumDistance(),
    )
    combined = interpolate(
        interpolation_model,
        kpoints,
        (
            BandEnergy(),
            BerryCurvature(win["fermi_energy"]; formulation = WYSV06()),
        ),
    )
    @test all(
        isapprox(
                view(combined.berry_curvature, :, :, index), ref_Ω[index];
                atol = 5.0e-5,
            )
            for index in eachindex(kpoints)
    )

    band_resolved = interpolate(
        interpolation_model,
        kpoints,
        BerryCurvature(; formulation = WYSV06BandResolved()),
    )
    for index in eachindex(kpoints)
        occupation = Int.(
            view(combined.band_energy, :, index) .<= win["fermi_energy"]
        )
        occupied_sum = dropdims(
            sum(
                view(band_resolved.berry_curvature, :, :, :, index) .*
                    reshape(occupation, 1, 1, :);
                dims = 3,
            );
            dims = 3,
        )
        @test isapprox(
            occupied_sum, view(combined.berry_curvature, :, :, index); atol = 1.0e-8
        )
    end

    separate = interpolate(
        interpolation_model,
        kpoints,
        BerryCurvature(win["fermi_energy"]; formulation = WYSV06()),
    )
    @test separate.berry_curvature == combined.berry_curvature
    lvts = interpolate(
        interpolation_model,
        kpoints,
        BerryCurvature(win["fermi_energy"]; formulation = LVTS12()),
    )
    @test all(
        isapprox(
                view(lvts.berry_curvature, :, :, index), ref_Ω[index];
                atol = 5.0e-5,
            )
            for index in eachindex(kpoints)
    )
    default_formulation = interpolate(
        interpolation_model, kpoints, BerryCurvature(win["fermi_energy"])
    )
    @test default_formulation.berry_curvature == lvts.berry_curvature
    @test_throws ArgumentError BerryCurvature(
        win["fermi_energy"]; formulation = WYSV06BandResolved()
    )
    @test_throws ArgumentError BerryCurvature(; formulation = WYSV06())
end
