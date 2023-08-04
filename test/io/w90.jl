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
        dataset"Si2_valence/reference/MDRS/Si2_valence", recip_lattice
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
    hamiltonian, position = read_w90_tb(dataset"Si2_valence/reference/WS/Si2_valence")

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
        dataset"Si2_valence/reference/MDRS/Si2_valence"
    )

    R1 = [-4, 0, 2]
    @test hamiltonian.Rvectors[1] == R1
    H111 = 0.0002681719333333333 + 5.364264333333333e-10im
    @test hamiltonian[1][1, 1] ≈ H111
    P11_origin = ComplexF64[0.67881607 + 0.0im, -0.67881621 + 0.0im, -0.67881622 + 0.0im]
    @test position[0, 0, 0][1, 1] ≈ P11_origin
end

@testitem "write_w90_tb" begin
    using Wannier.Datasets
    hamiltonian, position = Wannier.read_w90_tb(
        dataset"Si2_valence/reference/MDRS/Si2_valence"
    )

    outdir = mktempdir(; cleanup=true)
    outprefix = joinpath(outdir, "Si2_valence")
    write_w90_tb(outprefix, hamiltonian, position)

    hamiltonian2, position2 = Wannier.read_w90_tb(outprefix)
    @test hamiltonian ≈ hamiltonian2
    @test position ≈ position2
end

@testitem "read_w90_tb_chk_spn" begin
    using Wannier.Datasets
    using Wannier: Vec3

    hamiltonian, position, spin = Wannier.read_w90_tb_chk_spn(
        dataset"Fe_soc/reference/MDRS/Fe";
        chk=dataset"Fe_soc/reference/Fe.chk",
        spn=dataset"Fe_soc/Fe.spn",
    )

    S11 = Vec3(
        1.8048491616373703e-5 + 1.0545727604863155e-5im,
        -1.253138568448426e-5 + 1.1915668762749169e-5im,
        9.304292711693335e-6 + 5.54202611423551e-6im,
    )
    S12 = Vec3(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)
    S21 = Vec3(
        -2.954804706784439e-5 - 0.00013263850361096123im,
        -2.0567952787489034e-5 - 0.00014331577820929988im,
        -3.609044851399626e-5 + 5.085617072095425e-5im,
    )
    S22 = Vec3(
        0.00011145641163480431 + 3.648975777452525e-5im,
        -6.373032231840242e-5 + 5.626636785716623e-5im,
        -0.00011634213421945409 - 2.2934938576792313e-8im,
    )
    ref_S = [[S11, S21] [S12, S22]]
    @test spin[1][1:2, 1:2] ≈ ref_S
end
