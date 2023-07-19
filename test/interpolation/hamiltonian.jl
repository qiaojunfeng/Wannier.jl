@testitem "Hamiltonian WS" begin
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/reference/Si2_valence.chk"
    )
    hamiltonian = HamiltonianRspace(model; MDRS=false)

    ref_band = read_w90_band(dataset"Si2_valence/reference/ws/Si2_valence")
    eigenvalues = interpolate(ref_band.kpoints, hamiltonian)
    @test all(isapprox.(eigenvalues, ref_band.eigenvalues; atol=2e-5))
end

@testitem "Hamiltonian MDRS" begin
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/reference/Si2_valence.chk.fmt"
    )
    hamiltonian = HamiltonianRspace(model)

    # ref_band = read_w90_band(dataset"Si2_valence/reference/MDRS/Si2_valence")
    ref_band = read_w90_band(
        expanduser("~/git/WannierDatasets/datasets/Si2_valence/reference/MDRS/Si2_valence")
    )
    eigenvalues = interpolate(ref_band.kpoints, hamiltonian)
    @test all(isapprox.(eigenvalues, ref_band.eigenvalues; atol=2e-5))
end
