@testitem "read_w90_ibz constructs a localization-ready model" begin
    using Wannier.Datasets
    using LinearAlgebra

    mktempdir() do directory
        prefix = joinpath(directory, "Si2")
        sources = Dict(
            "nnkp" => dataset"Si2_hse/outputs/Si2.nnkp",
            "isym" => dataset"Si2_hse/Si2.isym",
            "ieig" => dataset"Si2_hse/Si2.ieig",
            "immn" => dataset"Si2_hse/Si2.immn",
            "iamn" => dataset"Si2_hse/Si2.iamn",
            "win" => dataset"Si2_hse/Si2.win",
        )
        for (extension, source) in sources
            symlink(source, "$prefix.$extension")
        end

        symmetric_model = read_w90_ibz(prefix; frozen_band_indices = 1:4)
        model = symmetric_model.model
        constraint = symmetric_model.constraint
        @test n_kpoints(model) == constraint.nk_fbz
        @test Wannier.n_kpoints_ibz(symmetric_model) == constraint.nk_ibz
        @test all(model.frozen_bands[1:4, :])
        @test !any(model.frozen_bands[5:end, :])
        @test size(symmetric_model.overlaps_ibz, 4) == constraint.nk_ibz
        gauge_ibz = Wannier.extract_ibz_gauges(model.gauges, constraint)
        @test maximum(
            opnorm(gauge_ibz[:, :, k]' * gauge_ibz[:, :, k] - I)
                for k in axes(gauge_ibz, 3)
        ) < 1.0e-12
        @test Wannier.covariance_residual(gauge_ibz, constraint) < 1.0e-7

        @test_throws ArgumentError read_w90_ibz(
            prefix;
            frozen_band_indices = 1:4,
            frozen_bands = model.frozen_bands,
        )
        @test_throws DimensionMismatch read_w90_ibz(
            prefix; frozen_bands = falses(1, 1)
        )
    end
end
