module BenchDisentangle

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2")
    U = model.gauges
    frozen_bands = model.frozen_bands
    nb = Wannier.n_bands(model)
    nw = Wannier.n_wannier(model)

    SUITE["U_to_X_Y"] = @benchmarkable Wannier.U_to_X_Y($U, $frozen_bands)

    X, Y = Wannier.U_to_X_Y(U, frozen_bands)
    SUITE["X_Y_to_U"] = @benchmarkable Wannier.X_Y_to_U($X, $Y)

    SUITE["X_Y_to_XY"] = @benchmarkable Wannier.X_Y_to_XY($X, $Y)

    XY = Wannier.X_Y_to_XY(X, Y)
    SUITE["XY_to_X_Y"] = @benchmarkable Wannier.XY_to_X_Y($XY, $nb, $nw)

    # end-to-end XYLayout path — 10 iterations
    SUITE["localize"] = @benchmarkable localize($model, max_iter = 10)

    # fused fg! closure, one call
    # Note: _make_fg! is private (for internal use) but exposed here for
    # low-level micro-benchmarking the kernel. Public API routes through localize().
    prob = Wannier.Problem(Wannier.Variance(), model)
    fg!  = Wannier._make_fg!(prob)
    XYbuf = copy(XY)
    Gbuf  = similar(XYbuf)
    SUITE["fg!"] = @benchmarkable $fg!(1.0, $Gbuf, $XYbuf)

end  # module

BenchDisentangle.SUITE
