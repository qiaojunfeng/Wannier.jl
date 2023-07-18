module BenchMaxloc

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

model = load_dataset("Si2_valence")

# just run 10 iterations
SUITE["max_localize"] = @benchmarkable max_localize($model, max_iter=10)

end  # module

BenchMaxloc.SUITE
