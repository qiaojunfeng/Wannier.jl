module BenchInterpolationConstruction

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2_valence")

    SUITE["Wigner-Seitz"] =
        @benchmarkable InterpolationModel($model; real_space = WignerSeitz())
    SUITE["minimum distance"] =
        @benchmarkable InterpolationModel($model; real_space = MinimumDistance())

end  # module

BenchInterpolationConstruction.SUITE
