module BenchMaxloc

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2_valence")

    # just run 10 iterations (UGauge path — isolated bands)
    SUITE["localize"] = @benchmarkable localize($model, max_iter = 10)

    # fused fg! closure, one call
    prob = Wannier.Problem(Wannier.Variance(), model)
    fg!  = Wannier._make_optim_fg!(prob)
    U    = copy(model.gauges)
    G    = similar(U)
    SUITE["fg!"] = @benchmarkable $fg!(1.0, $G, $U)

end  # module

BenchMaxloc.SUITE
