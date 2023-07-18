module BenchParalleltransport

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

model = load_dataset("Si2_valence")

SUITE["parallel_transport"] = @benchmarkable parallel_transport($model)

end  # module

BenchParalleltransport.SUITE
