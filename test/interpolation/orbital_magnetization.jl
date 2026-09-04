@testitem "OrbitalMagnetization observable" begin
    using LinearAlgebra
    using DelimitedFiles
    using Wannier.Datasets

    model = read_w90_with_chk(
        dataset"Fe_soc_coarse/Fe", dataset"Fe_soc_coarse/outputs/Fe.chk"
    )
    uHu = read_uHu(dataset"Fe_soc_coarse/Fe.uHu").uHu
    win = read_win(dataset"Fe_soc_coarse/Fe.win")
    ref_kpt = read_w90_band_kpt(
        dataset"Fe_soc_coarse/outputs/postw90/Fe-path.kpt"
    )
    ref_dat = readdlm(dataset"Fe_soc_coarse/outputs/postw90/Fe-morb.dat")
    # Wannier90 writes -1/2 of the LVTS12 Eq. 97 tensor.
    reference = map(eachrow(ref_dat[:, 2:end])) do value
        -2 * Wannier.axialvector_to_antisymmetrictensor(value)
    end

    segment = KSegment(
        reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"]
    )
    kpoints = collect(KPath(segment, 5))
    # postw90.x omits the H point at index 22.
    deleteat!(kpoints, 22)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1.0e-6)

    interpolation_model = InterpolationModel(
        model;
        operators = (;
            berry_connection = BerryConnection(; imlog_diag = false),
            hamiltonian_position = HamiltonianPosition(),
            position_hamiltonian_position = PositionHamiltonianPosition(uHu),
        ),
        real_space = MinimumDistance(),
    )
    result = interpolate(
        interpolation_model,
        kpoints,
        OrbitalMagnetization(win["fermi_energy"]),
    )
    @test all(
        isapprox(
                view(result.orbital_magnetization, :, :, index), reference[index];
                atol = 5.0e-7,
            )
            for index in eachindex(kpoints)
    )

    combined = interpolate(
        interpolation_model,
        kpoints,
        (
            BandEnergy(),
            BerryCurvature(win["fermi_energy"]),
            OrbitalMagnetization(win["fermi_energy"]),
        ),
    )
    @test combined.orbital_magnetization == result.orbital_magnetization
end
