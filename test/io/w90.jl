@testitem "read nnkp" begin
    using Wannier.Datasets
    kstencil = read_nnkp_compute_bweights(dataset"Si2_valence/reference/Si2_valence.nnkp")
    nnkp = read_nnkp(dataset"Si2_valence/reference/Si2_valence.nnkp.toml")

    @test reciprocal_lattice(kstencil) ≈ nnkp.recip_lattice
    @test kstencil.kpoints ≈ nnkp.kpoints
    @test kstencil.kpb_k == nnkp.kpb_k
    @test kstencil.kpb_G == nnkp.kpb_G

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
        3.361532, 3.361532, 3.361532, 3.361532, 3.361532, 3.361532, 3.361532, 3.361532
    ]
    @test isapprox(kstencil.bvectors, ref_bvectors; atol=1e-5)
    @test isapprox(kstencil.bweights, ref_bweights; atol=1e-5)
end

@testitem "read/write nnkp" begin
    using Wannier.Datasets
    kstencil = read_nnkp_compute_bweights(dataset"Si2_valence/reference/Si2_valence.nnkp")
    tmpfile = tempname(; cleanup=true)
    n_wann = 4
    write_nnkp(tmpfile, kstencil; n_wann)

    kstencil2 = read_nnkp_compute_bweights(tmpfile)
    @test kstencil ≈ kstencil2
end

@testitem "read/write w90 band" begin
    using Wannier.Datasets
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = reciprocal_lattice(win.unit_cell_cart)
    kpi, eigenvalues = read_w90_band(
        dataset"Si2_valence/reference/mdrs/Si2_valence", recip_lattice
    )

    outdir = mktempdir(; cleanup=true)
    outprefix = joinpath(outdir, "Si2_valence")
    write_w90_band(outprefix, kpi, eigenvalues)
    kpi2, eigenvalues2 = read_w90_band(outprefix, recip_lattice)

    @test kpi.kpaths ≈ kpi2.kpaths
    @test kpi.labels == kpi2.labels
    @test kpi.basis ≈ kpi2.basis
    @test Symbol(kpi.setting) == Symbol(kpi2.setting)
    @test eigenvalues ≈ eigenvalues2
end

@testitem "read_w90_tb WS" begin
    using Wannier.Datasets
    hamiltonian, position = read_w90_tb(dataset"Si2_valence/reference/ws/Si2_valence")

    # just some simple tests
    R1 = [-4, 0, 2]
    @test hamiltonian.Rvectors[1] == R1
    n_degen_R1 = 3
    H111 = 0.80451304E-03 + im * 0.16092791E-08
    @test hamiltonian[1][1, 1] ≈ H111 / n_degen_R1

    Rend = [4, 0, -2]
    @test position.Rvectors[end] == Rend
    P111end = [0.67881607 + 0.0im, -0.67881621 + 0.0im, -0.67881622 + 0.0im]
    @test position[0, 0, 0][1, 1] ≈ P111end
end

@testitem "read_w90_tb MDRS" begin
    using Wannier.Datasets
    hamiltonian, position = Wannier.read_w90_tb(
        dataset"Si2_valence/reference/mdrs/Si2_valence"
    )

    R1 = [-4, 0, 2]
    @test hamiltonian.Rvectors[1] == R1
    H111 = 0.0002681719333333333 + 5.364264333333333e-10im
    @test hamiltonian[1][1, 1] ≈ H111
    P11_origin = ComplexF64[0.67881607 + 0.0im, -0.67881621 + 0.0im, -0.67881622 + 0.0im]
    @test position[0, 0, 0][1, 1] ≈ P11_origin
end

@testitem "read_w90_tb_chk_spn" begin
    using Wannier.Datasets
    using Wannier: Vec3
    # hamiltonian, position, spin = Wannier.read_w90_tb_chk_spn(
    #     dataset"Fe/reference/Fe";
    #     spn=dataset"Fe/Fe.spn",
    # )
    hamiltonian, position, spin = Wannier.read_w90_tb_chk_spn(
        expanduser("~/git/WannierDatasets/datasets/Fe/reference/MDRS/Fe");
        chk=expanduser("~/git/WannierDatasets/datasets/Fe/reference/Fe.chk"),
        spn=expanduser("~/git/WannierDatasets/datasets/Fe/Fe.spn"),
    )

    S11 = Vec3(
        -4.648028910440591e-5 + 2.374737696435062e-5im,
        -6.0301956322061204e-5 - 1.638901806984331e-5im,
        -0.0001239322774898806 - 4.426989554163255e-7im,
    )
    S12 = Vec3(
        1.5735137180138543e-5 - 3.721098851980894e-5im,
        2.6168818950736308e-5 - 6.011199334516285e-5im,
        -3.1928160018176935e-5 + 8.37424892167478e-5im,
    )
    S21 = Vec3(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)
    S22 = Vec3(
        6.305695608742687e-6 + 1.0805021589902928e-6im,
        9.374556752615726e-6 - 4.398377968713008e-7im,
        -1.01162536765518e-5 - 1.7493947814005437e-7im,
    )
    ref_S = [[S11, S21] [S12, S22]]
    @test spin[1][1:2, 1:2] ≈ ref_S
end
