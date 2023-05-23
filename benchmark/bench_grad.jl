module BenchGrad

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))
bvectors = model.bvectors
M = model.M
U = model.U

SUITE["omega_grad"] = @benchmarkable omega_grad($bvectors, $M, $U)

end  # module

BenchGrad.SUITE
