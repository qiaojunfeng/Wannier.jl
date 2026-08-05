module BenchGrad

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2")
    kstencil = model.kstencil
    overlaps = model.overlaps
    gauges = model.gauges

    ws = Wannier.Workspace(model)
    Wannier.compute_MU_UtMU!(ws, kstencil, overlaps, gauges)
    SUITE["omega_grad!"] = @benchmarkable Wannier.omega_grad!($(ws.GU), $ws, $kstencil, $overlaps)

end  # module

BenchGrad.SUITE
