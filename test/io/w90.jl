@testitem "read nnkp" begin
    using Wannier.Datasets
    kstencil = read_nnkp_compute_bweights(dataset"Si2_valence/outputs/Si2_valence.nnkp")
    nnkp = read_nnkp(dataset"Si2_valence/outputs/Si2_valence.nnkp.toml")

    @test reciprocal_lattice(kstencil) ≈ nnkp["recip_lattice"]
    @test kstencil.kpoints ≈ nnkp["kpoints"]
    @test kstencil.kpb_k == nnkp["kpb_k"]
    @test kstencil.kpb_G == nnkp["kpb_G"]

    # copied from wout file
    ref_bvectors = [
        [0.192835, 0.192835, -0.192835],
        [0.192835, -0.192835, 0.192835],
        [-0.192835, 0.192835, 0.192835],
        [0.192835, 0.192835, 0.192835],
        [-0.192835, -0.192835, 0.192835],
        [-0.192835, 0.192835, -0.192835],
        [0.192835, -0.192835, -0.192835],
        [-0.192835, -0.192835, -0.192835],
    ]
    ref_bweights = [
        3.361532, 3.361532, 3.361532, 3.361532, 3.361532, 3.361532, 3.361532, 3.361532,
    ]
    @test isapprox(kstencil.bvectors, ref_bvectors; atol = 1.0e-5)
    @test isapprox(kstencil.bweights, ref_bweights; atol = 1.0e-5)
end

@testitem "read/write nnkp" begin
    using Wannier.Datasets
    kstencil = read_nnkp_compute_bweights(dataset"Si2_valence/outputs/Si2_valence.nnkp")
    tmpfile = tempname(; cleanup = true)
    n_wann = 4
    write_nnkp(tmpfile, kstencil; n_wann)

    kstencil2 = read_nnkp_compute_bweights(tmpfile)
    @test kstencil ≈ kstencil2
end

@testitem "read_w90_tb WS" begin
    using Wannier.Datasets
    model = read_w90_tb(dataset"Si2_valence/outputs/WS/Si2_valence")

    # just some simple tests
    R1 = Vec3(-4, 0, 2)
    R1_index = model.real_space.vector_index[R1]
    n_degen_R1 = 3
    H111 = 0.80451304e-3 + im * 0.16092791e-8
    @test model.operators.hamiltonian.coefficients[1, 1, R1_index] ≈
        H111 / n_degen_R1

    Rend = Vec3(4, 0, -2)
    @test haskey(model.real_space.vector_index, Rend)
    P111end = [0.67881607 + 0.0im, -0.67881621 + 0.0im, -0.67881622 + 0.0im]
    origin_index = model.real_space.vector_index[Vec3(0, 0, 0)]
    @test model.operators.berry_connection.coefficients[1, 1, :, origin_index] ≈
        P111end
end

@testitem "read_w90_tb MDRS" begin
    using Wannier.Datasets
    model = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")

    R1 = Vec3(-4, 0, 2)
    R1_index = model.real_space.vector_index[R1]
    H111 = 0.0002681719333333333 + 5.364264333333333e-10im
    @test model.operators.hamiltonian.coefficients[1, 1, R1_index] ≈ H111
    P11_origin = ComplexF64[0.67881607 + 0.0im, -0.67881621 + 0.0im, -0.67881622 + 0.0im]
    origin_index = model.real_space.vector_index[Vec3(0, 0, 0)]
    @test model.operators.berry_connection.coefficients[1, 1, :, origin_index] ≈
        P11_origin
end

@testitem "write_w90_tb" begin
    using Wannier.Datasets
    model = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")

    outdir = mktempdir(; cleanup = true)
    outprefix = joinpath(outdir, "Si2_valence")
    write_w90_tb(outprefix, model)

    model2 = read_w90_tb(outprefix)
    @test model.real_space.vectors == model2.real_space.vectors
    @test model.basis.fractional_centers ≈ model2.basis.fractional_centers
    @test model.operators.hamiltonian.coefficients ≈
        model2.operators.hamiltonian.coefficients
    @test model.operators.berry_connection.coefficients ≈
        model2.operators.berry_connection.coefficients
end

@testitem "read_w90_hr WS" begin
    using Wannier.Datasets
    reference = read_w90_tb(dataset"Si2_valence/outputs/WS/Si2_valence")
    model = read_w90_hr(
        dataset"Si2_valence/outputs/WS/Si2_valence",
        real_lattice(reference);
        fractional_centers = reference.basis.fractional_centers,
    )
    kpoints = [Vec3(0.13, 0.21, -0.07), Vec3(-0.19, 0.04, 0.31)]
    energy = interpolate(model, kpoints, BandEnergy()).band_energy
    reference_energy = interpolate(
        reference, kpoints, BandEnergy()
    ).band_energy
    @test maximum(abs, energy - reference_energy) < 5.0e-5
end

@testitem "read_w90_hr MDRS" begin
    using Wannier.Datasets
    reference = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")
    model = read_w90_hr(
        dataset"Si2_valence/outputs/MDRS/Si2_valence",
        real_lattice(reference);
        fractional_centers = reference.basis.fractional_centers,
    )
    kpoints = [Vec3(0.13, 0.21, -0.07), Vec3(-0.19, 0.04, 0.31)]
    energy = interpolate(model, kpoints, BandEnergy()).band_energy
    reference_energy = interpolate(
        reference, kpoints, BandEnergy()
    ).band_energy
    @test maximum(abs, energy - reference_energy) < 5.0e-5
end

@testitem "write_w90_hr" begin
    using Wannier.Datasets
    reference = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")

    outdir = mktempdir(; cleanup = true)
    outprefix = joinpath(outdir, "Si2_valence")
    write_w90_hr(outprefix, reference)

    model = read_w90_hr(
        outprefix,
        real_lattice(reference);
        fractional_centers = reference.basis.fractional_centers,
    )
    @test model.real_space.vectors == reference.real_space.vectors
    @test maximum(
        abs,
        reference.operators.hamiltonian.coefficients -
            model.operators.hamiltonian.coefficients,
    ) < 1.0e-6
end

@testitem "read_w90_tb_chk_spn" begin
    using Wannier.Datasets

    model = Wannier.read_w90_tb_chk_spn(
        dataset"Fe_soc_coarse/outputs/Fe";
        chk = dataset"Fe_soc_coarse/outputs/Fe.chk",
        spn = dataset"Fe_soc_coarse/Fe.spn",
    )

    source_model = read_w90_with_chk(
        dataset"Fe_soc_coarse/Fe",
        dataset"Fe_soc_coarse/outputs/Fe.chk",
    )
    reference = InterpolationModel(
        source_model;
        operators = (;
            spin = BlochOperator(read_spn(dataset"Fe_soc_coarse/Fe.spn")),
        ),
        real_space = MinimumDistance(),
    )
    kpoints = [Vec3(0.13, 0.21, -0.07), Vec3(-0.19, 0.04, 0.31)]
    spin = interpolate(model, kpoints, SpinExpectation()).spin_expectation
    reference_spin = interpolate(
        reference, kpoints, SpinExpectation()
    ).spin_expectation
    @test maximum(abs, spin - reference_spin) < 5.0e-7
end
