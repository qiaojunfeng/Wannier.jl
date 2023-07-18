module BenchBvector

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

win = read_win(dataset"Si2/Si2.win")
kpoints = win.kpoints
recip_lattice = get_recip_lattice(win.unit_cell_cart)

SUITE["get_bvectors"] = @benchmarkable get_bvectors($kpoints, $recip_lattice)

end  # module

BenchBvector.SUITE
