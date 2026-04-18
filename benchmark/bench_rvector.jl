module BenchRvectors

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    win = read_win(dataset"Si2_valence/Si2_valence.win")
    lattice = win["unit_cell_cart"]
    kgrid = win["mp_grid"]

    wout = read_wout(dataset"Si2_valence/outputs/Si2_valence.wout")
    # to fractional coordinates for MDRS T-vector search
    centers = map(c -> inv(lattice) * c, wout["centers"])

    SUITE["WignerSeitzRspace"] = @benchmarkable Wannier.WignerSeitzRspace($lattice, $kgrid)
    SUITE["MDRSRspace"] =
        @benchmarkable Wannier.MDRSRspace($lattice, $kgrid, $centers)

end  # module

BenchRvectors.SUITE
