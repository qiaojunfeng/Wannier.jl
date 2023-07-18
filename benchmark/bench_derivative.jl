module BenchDerivative

using BenchmarkTools
using Wannier
using Wannier: Vec3
using Wannier.Datasets

SUITE = BenchmarkGroup()

tb = read_w90_tb(dataset"Si2_valence/reference/mdrs/Si2_valence")
Rvectors = tb.Rvectors
H = tb.H
# choose a random k-point such that there is no degeneracy
k = [Vec3(0.1, 0.2, 0.3)]

SUITE["velocity_fd"] = @benchmarkable Wannier.velocity_fd(Rvectors, H, k)
SUITE["velocity"] = @benchmarkable Wannier.velocity(H, k)
SUITE["get_dH_da"] = @benchmarkable Wannier.get_dH_da(H, k)
SUITE["effmass_fd"] = @benchmarkable Wannier.effmass_fd(Rvectors, H, k)
SUITE["get_d2H_dadb"] = @benchmarkable Wannier.get_d2H_dadb(H, k)

end  # module

BenchDerivative.SUITE
