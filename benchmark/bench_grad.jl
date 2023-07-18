module BenchGrad

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

model = load_dataset("Si2")
bvectors = model.bvectors
M = model.M
U = model.U

SUITE["omega_grad"] = @benchmarkable omega_grad($bvectors, $M, $U)

end  # module

BenchGrad.SUITE
