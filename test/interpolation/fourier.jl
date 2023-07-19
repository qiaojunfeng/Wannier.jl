@testitem "fourier" begin
    using Wannier.Datasets
    tbdat = read_w90_tbdat(dataset"Si2_valence/reference/ws/Si2_valence_tb.dat")
    Rdomain = Wannier.WSRspaceDomain(tbdat.lattice, tbdat.Rvectors, tbdat.Rdegens)
    ref_Hᴿ = tbdat.H

    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/reference/Si2_valence.chk.fmt"
    )
    Hᵏ = HamiltonianKspace(model)
    Hᴿ = fourier(Hᵏ, Rdomain)
    @test all(isapprox.(ref_Hᴿ, Hᴿ; atol=5e-7))
end

@testitem "invfourier" begin
    using Wannier.Datasets
    hamiltonian, position = read_w90_tb(dataset"Si2_valence/reference/ws/Si2_valence")
    ref_band = read_w90_band(dataset"Si2_valence/reference/ws/Si2_valence")
    klist = Wannier.KpointList(
        reciprocal_lattice(real_lattice(hamiltonian)), ref_band.kpoints
    )

    Hᵏ = invfourier(hamiltonian, klist)
    hamiltonian_k = HamiltonianKspace(klist, Hᵏ)
    eigenvals, eigenvecs = eigen(hamiltonian_k)
    @test all(isapprox.(eigenvals, ref_band.eigenvalues; atol=2e-5))
end
