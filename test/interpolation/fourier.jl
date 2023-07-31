@testitem "fourier" begin
    using Wannier.Datasets
    tbdat = read_w90_tbdat(dataset"Si2_valence/reference/WS/Si2_valence_tb.dat")
    ref_Hᴿ = tbdat.H
    Rspace = Wannier.WignerSeitzRspace(tbdat.lattice, tbdat.Rvectors, tbdat.Rdegens)

    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/reference/Si2_valence.chk.fmt"
    )
    Hᵏ = transform_gauge(model.eigenvalues, model.gauges)
    Hᴿ = fourier(model.kpoints, Hᵏ, Rspace)
    @test all(isapprox.(ref_Hᴿ, Hᴿ; atol=5e-7))
end

@testitem "invfourier" begin
    using Wannier.Datasets
    hamiltonian, position = read_w90_tb(dataset"Si2_valence/reference/WS/Si2_valence")
    ref_band = read_w90_band(dataset"Si2_valence/reference/WS/Si2_valence")
    kpoints = ref_band.kpoints

    Hᵏ = invfourier(hamiltonian, kpoints)
    eigenvals, eigenvecs = eigen(Hᵏ)
    @test all(isapprox.(eigenvals, ref_band.eigenvalues; atol=2e-5))
end
