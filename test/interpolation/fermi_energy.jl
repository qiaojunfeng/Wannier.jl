@testitem "compute_fermi_energy" begin
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_valence/Si2_valence", dataset"Si2_valence/reference/Si2_valence.chk.fmt"
    )
    hamiltonian = TBHamiltonian(model)
    interp = HamiltonianInterpolator(hamiltonian)
    kgrid = [12, 12, 12]

    εF = Wannier.compute_fermi_energy(kgrid, interp, 7.1, 1e-2, Wannier.ColdSmearing())
    εF_ref = 4.637512665199997
    @test isapprox(εF, εF_ref; atol=1e-3)
end

@testitem "compute_fermi_energy graphene" begin
    using Wannier.Datasets
    model = load_dataset("graphene")
    model.gauges .= read_amn(dataset"graphene/reference/graphene.dis.amn")
    hamiltonian = TBHamiltonian(model)
    interp = HamiltonianInterpolator(hamiltonian)
    # on purposely choose 5x5x1 since this grid skips the K point, and
    # a simple Fermi energy search would fail.
    kgrid = [5, 5, 1]

    εF = Wannier.compute_fermi_energy(kgrid, interp, 8, 0, Wannier.NoneSmearing())
    εF_ref = -1.03673405699654
    @test isapprox(εF, εF_ref; atol=1e-3)
end
