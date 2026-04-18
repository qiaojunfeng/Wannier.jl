module BenchBvector

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    win = read_win(dataset"Si2/Si2.win")
    kpoints = win["kpoints"]
    mp_grid = win["mp_grid"]
    recip_lattice = reciprocal_lattice(win["unit_cell_cart"])

    SUITE["generate_kspace_stencil"] =
        @benchmarkable generate_kspace_stencil($recip_lattice, $mp_grid, $kpoints)

end  # module

BenchBvector.SUITE
