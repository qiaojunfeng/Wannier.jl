module BenchDerivative

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))
Rvecs = model.kRvectors.Rvectors
H = model.H
# choose a random k-point such that there is no degeneracy
k = reshape([0.1, 0.2, 0.3], (3, 1))

SUITE["velocity_fd"] = @benchmarkable Wannier.velocity_fd(Rvecs, H, k)
SUITE["velocity"] = @benchmarkable Wannier.velocity(Rvecs, H, k)
SUITE["get_dH_da"] = @benchmarkable Wannier.get_dH_da(Rvecs, H, k)
SUITE["effmass_fd"] = @benchmarkable Wannier.effmass_fd(Rvecs, H, k)
SUITE["get_d2H_dadb"] = @benchmarkable Wannier.get_d2H_dadb(Rvecs, H, k)

end  # module

BenchDerivative.SUITE
