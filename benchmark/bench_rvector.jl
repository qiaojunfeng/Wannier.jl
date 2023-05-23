module BenchRvectors

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
lattice = win.unit_cell_cart
mp_grid = win.mp_grid

wout = read_wout(joinpath(FIXTURE_PATH, "valence/band/silicon.wout"))
# to fractional coordinates
centers = inv(lattice) * wout.centers

SUITE["get_Rvectors_ws"] = @benchmarkable Wannier.get_Rvectors_ws($lattice, $mp_grid)
SUITE["get_Rvectors_mdrs"] = @benchmarkable Wannier.get_Rvectors_mdrs(
    $lattice, $mp_grid, $centers
)

end  # module

BenchRvectors.SUITE
