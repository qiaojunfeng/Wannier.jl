module BenchRvectors

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

win = read_win(dataset"Si2_valence/Si2_valence.win")
lattice = win.unit_cell_cart
kgrid = win.mp_grid

wout = read_wout(dataset"Si2_valence/reference/Si2_valence.wout")
# to fractional coordinates
centers = map(c -> inv(lattice) * c, wout.centers)

SUITE["get_Rvectors_ws"] = @benchmarkable Wannier.get_Rvectors_ws($lattice, $kgrid)
SUITE["get_Rvectors_mdrs"] = @benchmarkable Wannier.get_Rvectors_mdrs(
    $lattice, $kgrid, $centers
)

end  # module

BenchRvectors.SUITE
