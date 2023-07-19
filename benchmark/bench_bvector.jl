module BenchBvector

using BenchmarkTools
using Wannier
using Wannier.Datasets

SUITE = BenchmarkGroup()

win = read_win(dataset"Si2/Si2.win")
kpoints = win.kpoints
recip_lattice = reciprocal_lattice(win.unit_cell_cart)

SUITE["generate_bvectors"] = @benchmarkable generate_bvectors($kpoints, $recip_lattice)

end  # module

BenchBvector.SUITE
