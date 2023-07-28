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
