module BenchSpread

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))
bvectors = model.bvectors
M = model.M
U = model.U

SUITE["omega"] = @benchmarkable omega($bvectors, $M, $U)

end  # module

BenchSpread.SUITE
