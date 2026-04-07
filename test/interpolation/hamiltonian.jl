@testitem "Hamiltonian WS" begin
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/outputs/Si2_valence.chk.fmt"
    )
    hamiltonian = TBHamiltonian(model; MDRS = false)
    interp = HamiltonianInterpolator(hamiltonian)

    ref_band = read_w90_band(dataset"Si2_valence/outputs/WS/Si2_valence")
    # if I use the kpoints in ref_band, the difference between eigenvalues is
    # around 2e-5, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_band.kpoints
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    kseg = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    kpath = KPath(kseg)
    eigenvalues = interp(collect(kpath))[1]
    @test all(isapprox.(eigenvalues, ref_band.eigenvalues; atol = 2.0e-6))
end

@testitem "Hamiltonian MDRS" begin
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/outputs/Si2_valence.chk.fmt"
    )
    hamiltonian = TBHamiltonian(model)
    interp = HamiltonianInterpolator(hamiltonian)

    ref_band = read_w90_band(dataset"Si2_valence/outputs/MDRS/Si2_valence")

    # if I use the kpoints in ref_band, the difference between eigenvalues is
    # around 2e-5, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_band.kpoints
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    kseg = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    kpath = KPath(kseg)
    eigenvalues = interp(collect(kpath))[1]
    @test all(isapprox.(eigenvalues, ref_band.eigenvalues; atol = 1.0e-7))
end
