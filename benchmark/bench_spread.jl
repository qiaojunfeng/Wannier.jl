module BenchSpread

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2")
    kstencil = model.kstencil
    overlaps = model.overlaps
    gauges = model.gauges

    SUITE["spread"] = @benchmarkable spread($kstencil, $overlaps, $gauges)

    ws = Wannier.Workspace(model)
    Wannier.compute_MU_UtMU!(ws, kstencil, overlaps, gauges)
    SUITE["omega!"] = @benchmarkable Wannier.omega!($ws, $kstencil, $overlaps)

end  # module

BenchSpread.SUITE
