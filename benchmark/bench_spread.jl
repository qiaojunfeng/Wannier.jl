module BenchSpread

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

model = load_dataset("Si2")
bvectors = model.bvectors
M = model.M
U = model.U

SUITE["omega"] = @benchmarkable omega($bvectors, $M, $U)

end  # module

BenchSpread.SUITE
