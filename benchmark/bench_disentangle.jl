module BenchDisentangle

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))
U = model.U
frozen_bands = model.frozen_bands

SUITE["U_to_X_Y"] = @benchmarkable Wannier.U_to_X_Y($U, $frozen_bands)

X, Y = Wannier.U_to_X_Y(model.U, model.frozen_bands)
SUITE["X_Y_to_U"] = @benchmarkable Wannier.X_Y_to_U($X, $Y)

SUITE["XY_to_X_Y"] = @benchmarkable Wannier.X_Y_to_XY($X, $Y)

XY = Wannier.X_Y_to_XY(X, Y)
n_bands = model.n_bands
n_wann = model.n_wann
SUITE["X_Y_to_XY"] = @benchmarkable Wannier.XY_to_X_Y($XY, $n_bands, $n_wann)

# just run 10 iterations
SUITE["disentangle"] = @benchmarkable disentangle(SpreadPenalty(), $model, max_iter=10)

end  # module

BenchDisentangle.SUITE
