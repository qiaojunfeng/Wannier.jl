module BenchParalleltransport

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90(joinpath(FIXTURE_PATH, "valence", "silicon"))

SUITE["parallel_transport"] = @benchmarkable parallel_transport($model)

end  # module

BenchParalleltransport.SUITE
