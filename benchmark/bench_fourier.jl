module BenchFourier

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

model = load_dataset("Si2_valence")
model.U .= get_U(read_chk(dataset"Si2_valence/reference/Si2_valence.chk.fmt"))
tb_ws = read_w90_tb(dataset"Si2_valence/reference/ws/Si2_valence")
tb_mdrs = read_w90_tb(dataset"Si2_valence/reference/mdrs/Si2_valence")

Hᵏ = Wannier.get_Hk(model.E, model.U)
kRvectors_ws = tb_ws.Rvectors
kRvectors_mdrs = tb_mdrs.Rvectors

SUITE["get_Hk"] = @benchmarkable Wannier.get_Hk($(model.E), $(model.U))

tbhami = tb_ws.H

SUITE["fourier WS"] = @benchmarkable Wannier.HR_ws(
    $Hᵏ, $(model.kpoints), $(kRvectors_ws.R), $(length(model.E[1]))
)

SUITE["fourier MDRS v2"] = @benchmarkable (
    HR = Wannier.HR_ws($Hᵏ, $(model.kpoints), $(kRvectors_ws.R), $(length(model.E[1])));
    Wannier.mdrs_v1tov2(HR, $(kRvectors_mdrs))
)

out = [Wannier.zeros_block(tbhami) for i in 1:length(model.kpoints)]
SUITE["invfourier MDRS v2"] = @benchmarkable map(enumerate($(model.kpoints))) do (i, k)
    Wannier.invfourier($tbhami, k) do ib, iR, R_cart, b, fac
        @inbounds $(out)[i][ib] += fac * b.block[ib]
    end
end

end  # module

BenchFourier.SUITE
