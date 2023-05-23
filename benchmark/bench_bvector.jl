module BenchBvector

import Main.FIXTURE_PATH
using BenchmarkTools
using Wannier

SUITE = BenchmarkGroup()

win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))
kpoints = win.kpoints
recip_lattice = get_recip_lattice(win.unit_cell_cart)

SUITE["get_bvectors"] = @benchmarkable get_bvectors($kpoints, $recip_lattice)

end  # module

BenchBvector.SUITE
