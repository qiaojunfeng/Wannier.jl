module BenchFourier

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    # read_w90_tb already calls simplify -> BareRspace for fast invfourier
    ham_ws, _   = read_w90_tb(dataset"Si2_valence/outputs/WS/Si2_valence")
    ham_mdrs, _ = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")

    model = load_dataset("Si2_valence")
    kpoints = Wannier.kpoints(model)

    SUITE["invfourier WS"]   = @benchmarkable invfourier($ham_ws, $kpoints)
    SUITE["invfourier MDRS"] = @benchmarkable invfourier($ham_mdrs, $kpoints)

    # in-place variants reuse a preallocated output buffer
    Hk_ws_buf   = [zero(ham_ws[1])   for _ in 1:length(kpoints)]
    Hk_mdrs_buf = [zero(ham_mdrs[1]) for _ in 1:length(kpoints)]
    SUITE["invfourier! WS"]   = @benchmarkable Wannier.invfourier!($Hk_ws_buf,   $ham_ws,   $kpoints)
    SUITE["invfourier! MDRS"] = @benchmarkable Wannier.invfourier!($Hk_mdrs_buf, $ham_mdrs, $kpoints)

end  # module

BenchFourier.SUITE
