# A reusable fixture for a model
# no disentanglement

@testitem "exchanges" begin
    using LinearAlgebra
    using Wannier: Vec3
    using Wannier.Datasets
    model_up = read_w90_with_chk(dataset"Fe_collinear/Fe_up", dataset"Fe_collinear/reference/Fe_up.chk")
    model_dn = read_w90_with_chk(dataset"Fe_collinear/Fe_dn", dataset"Fe_collinear/reference/Fe_dn.chk")
    Mupdn = read_amn(dataset"Fe_collinear/Fe_updn.mud")
    model = Wannier.MagModel(model_up, model_dn, Mupdn)
    hami = Wannier.TBHamiltonian(model)

    atoms = Wannier.Atoms(model)

    n_ωh = 100
    nk = (5, 5, 5)
    exch = Wannier.calc_exchanges(hami, atoms, model_up.lattice, 0.01; R = Vec3(0,0,1), n_ωh=n_ωh, nk=nk)

    @test isapprox(sum(exch[1].J), 0.02541964326559294)
end
