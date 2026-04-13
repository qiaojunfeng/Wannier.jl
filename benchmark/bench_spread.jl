module BenchSpread

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2")
    kstencil = model.kstencil
    overlaps = model.overlaps
    gauges = model.gauges

    SUITE["omega"] = @benchmarkable omega($kstencil, $overlaps, $gauges)

end  # module

BenchSpread.SUITE
