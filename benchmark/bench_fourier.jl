module BenchFourier

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90(joinpath(FIXTURE_PATH, "valence/band/silicon"); amn=false)
model.U .= get_U(read_chk(joinpath(FIXTURE_PATH, "valence/band/silicon.chk.fmt")))
model_ws = read_w90_tb(
    joinpath(FIXTURE_PATH, "valence/band/ws/silicon"); kpoints=model.kpoints
)
model_mdrs = read_w90_tb(
    joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"); kpoints=model.kpoints
)

Hᵏ = Wannier.get_Hk(model.E, model.U)
kRvectors_ws = model_ws.kRvectors
kRvectors_mdrs = model_mdrs.kRvectors

SUITE["get_Hk"] = @benchmarkable Wannier.get_Hk($(model.E), $(model.U))

SUITE["fourier WS"] = @benchmarkable Wannier.fourier($kRvectors_ws, $Hᵏ)

SUITE["invfourier WS"] = @benchmarkable Wannier.invfourier(
    $kRvectors_ws, $(model_ws.H), $(model.kpoints)
)

SUITE["fourier MDRS v1"] = @benchmarkable Wannier.fourier($kRvectors_mdrs, $Hᵏ; version=:v1)

SUITE["invfourier MDRS v1"] = @benchmarkable Wannier.invfourier(
    $kRvectors_mdrs, $(model_ws.H), $(model.kpoints); version=:v1
)

SUITE["fourier MDRS v2"] = @benchmarkable Wannier.fourier($kRvectors_mdrs, $Hᵏ; version=:v2)

Hᴿ = Wannier.fourier(kRvectors_mdrs, Hᵏ; version=:v2)
SUITE["invfourier MDRS v2"] = @benchmarkable Wannier.invfourier(
    $kRvectors_mdrs, $Hᴿ, $(model.kpoints); version=:v2
)

end  # module

BenchFourier.SUITE
