module BenchDisentangle

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2")
    layout = Wannier.XYLayout()
    SUITE["initial XY parameters"] =
        @benchmarkable Wannier.initial_parameters($layout, $model)

    XY = Wannier.initial_parameters(layout, model)
    workspace = Wannier.Workspace(model)
    SUITE["assemble compact XY gauge"] =
        @benchmarkable Wannier.assemble_gauge!($layout, $XY, $model, $workspace)
    Wannier.assemble_gauge!(layout, XY, model, workspace)
    workspace.GU .= workspace.U
    G = similar(XY)
    SUITE["pull back compact XY gradient"] =
        @benchmarkable Wannier.pullback_gradient!($G, $layout, $model, $workspace)

    # end-to-end XYLayout path — 10 iterations
    SUITE["localize"] = @benchmarkable localize($model, max_iter = 10)

    # fused fg! closure, one call
    # Note: `_optimizer_callback` is private (for internal use) but exposed here for
    # low-level micro-benchmarking the kernel. Public API routes through localize().
    prob = Wannier.Problem(Wannier.Variance(), model)
    fg! = Wannier._optimizer_callback(prob)
    XYbuf = copy(XY)
    SUITE["fg!"] = @benchmarkable $fg!(1.0, $G, $XYbuf)

end  # module

BenchDisentangle.SUITE
