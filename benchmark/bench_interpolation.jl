module BenchInterpolation

    using BenchmarkTools
    using Wannier
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    model = load_dataset("Si2_valence")
    kpoints = Wannier.kpoints(model)

    wigner_seitz_model = InterpolationModel(model; real_space = WignerSeitz())
    minimum_distance_model = InterpolationModel(model; real_space = MinimumDistance())

    # Warm up the complete public path before BenchmarkTools measures it.
    interpolate(wigner_seitz_model, kpoints, BandEnergy())
    interpolate(minimum_distance_model, kpoints, BandEnergy())

    SUITE["band energy"]["Wigner-Seitz"] =
        @benchmarkable interpolate($wigner_seitz_model, $kpoints, BandEnergy())
    SUITE["band energy"]["minimum distance"] =
        @benchmarkable interpolate($minimum_distance_model, $kpoints, BandEnergy())

end  # module

BenchInterpolation.SUITE
