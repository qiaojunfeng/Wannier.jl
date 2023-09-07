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

    n_ωh = 300
    n_ωv = 50
    nk = (5, 5, 5)
    exch = Wannier.calc_exchanges(hami, atoms, model.lattice, 0.01; R = Vec3(0,0,1), n_ωh, n_ωv, nk)

    @code_warntype Wannier.ExchangeKGrid(hami, Wannier.uniform_shifted_kgrid(nk...), R)
    maxJ = maximum([tr(e.J) for e in exch])
    @test isapprox(maxJ, -158.90552565580805)
end
