module BenchFourier

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90(joinpath(FIXTURE_PATH, "valence/band/silicon"); amn=false)
model.U .= get_U(read_chk(joinpath(FIXTURE_PATH, "valence/band/silicon.chk.fmt")))
model_ws = read_w90_tb(
    joinpath(FIXTURE_PATH, "valence/band/ws/silicon")
)
model_mdrs = read_w90_tb(
    joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon")
)

Hᵏ = Wannier.get_Hk(model.E, model.U)
kRvectors_ws = model_ws.R
kRvectors_mdrs = model_mdrs.R

SUITE["get_Hk"] = @benchmarkable Wannier.get_Hk($(model.E), $(model.U))

tbhami = model_ws.H


SUITE["fourier WS"] = @benchmarkable Wannier.HR_ws($Hᵏ, $(model.kpoints), $(kRvectors_ws.R), $(length(model.E[1])))

SUITE["fourier MDRS v2"] = @benchmarkable (HR = Wannier.HR_ws($Hᵏ, $(model.kpoints), $(kRvectors_ws.R), $(length(model.E[1]))); Wannier.mdrs_v1tov2(HR, $(kRvectors_mdrs)))

out = [Wannier.zeros_block(tbhami) for i = 1:length(model.kpoints)]
SUITE["invfourier MDRS v2"] = @benchmarkable map(enumerate($(model.kpoints))) do (i, k)
    Wannier.invfourier($tbhami, k) do ib, iR, R_cart, b, fac
        @inbounds $(out)[i][ib] += fac * b.block[ib]
    end
end


end  # module

BenchFourier.SUITE
