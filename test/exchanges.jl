using LinearAlgebra
using Wannier: Vec3
# A reusable fixture for a model
# no disentanglement
model_up = read_w90(joinpath(FIXTURE_PATH, "iron", "up"))
model_dn = read_w90(joinpath(FIXTURE_PATH, "iron", "dn"))
Mupdn = read_amn(joinpath(FIXTURE_PATH, "iron", "iron.mud"))
model = Wannier.MagModel(model_up, model_dn, Mupdn)

@testset "exchanges" begin
    hami = Wannier.TBHamiltonian(model)

    atoms = Wannier.Atoms(model)

    n_ωh = 300
    n_ωv = 50
    nk = (5, 5, 5)
    exch = Wannier.calc_exchanges(hami, atoms, 0.01; R = Vec3(0,0,1), n_ωh, n_ωv, nk)

    @code_warntype Wannier.ExchangeKGrid(hami, Wannier.uniform_shifted_kgrid(nk...), R)
    maxJ = maximum([tr(e.J) for e in exch])
    @test isapprox(maxJ, -158.90552565580805)
end
