module BenchMaxloc

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90(joinpath(FIXTURE_PATH, "valence/silicon"))

# just run 10 iterations
SUITE["max_localize"] = @benchmarkable max_localize($model, max_iter=10)

end  # module

BenchMaxloc.SUITE
