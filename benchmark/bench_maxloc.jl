module BenchMaxloc

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2_valence")

    # just run 10 iterations (ULayout path — isolated bands)
    SUITE["localize"] = @benchmarkable localize($model, max_iter = 10)

    # fused fg! closure, one call
    # Note: `_optimizer_callback` is private (for internal use) but exposed here for
    # low-level micro-benchmarking the kernel. Public API routes through localize().
    prob = Wannier.Problem(Wannier.Variance(), model)
    fg! = Wannier._optimizer_callback(prob)
    U = copy(model.gauges)
    G = similar(U)
    SUITE["fg!"] = @benchmarkable $fg!(1.0, $G, $U)

end  # module

BenchMaxloc.SUITE
