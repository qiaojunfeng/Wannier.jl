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
    hamiltonian, position = read_w90_tb(dataset"Si2_valence/outputs/WS/Si2_valence")

    # just some simple tests
    R1 = [-4, 0, 2]
    @test hamiltonian.Rvectors[1] == R1
    n_degen_R1 = 3
    H111 = 0.80451304e-3 + im * 0.16092791e-8
    @test hamiltonian[1][1, 1] ≈ H111 / n_degen_R1

    Rend = [4, 0, -2]
    @test position.Rvectors[end] == Rend
    P111end = [0.67881607 + 0.0im, -0.67881621 + 0.0im, -0.67881622 + 0.0im]
    @test position[0, 0, 0][1, 1] ≈ P111end
end

@testitem "read_w90_tb MDRS" begin
    using Wannier.Datasets
    hamiltonian, position = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")

    R1 = [-4, 0, 2]
    @test hamiltonian.Rvectors[1] == R1
    H111 = 0.0002681719333333333 + 5.364264333333333e-10im
    @test hamiltonian[1][1, 1] ≈ H111
    P11_origin = ComplexF64[0.67881607 + 0.0im, -0.67881621 + 0.0im, -0.67881622 + 0.0im]
    @test position[0, 0, 0][1, 1] ≈ P11_origin
end

@testitem "write_w90_tb" begin
    using Wannier.Datasets
    hamiltonian, position = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")

    outdir = mktempdir(; cleanup = true)
    outprefix = joinpath(outdir, "Si2_valence")
    write_w90_tb(outprefix, hamiltonian, position)

    hamiltonian2, position2 = read_w90_tb(outprefix)
    @test hamiltonian ≈ hamiltonian2
    @test position ≈ position2
end

@testitem "read_w90_hr WS" begin
    using Wannier.Datasets
    ref_hamiltonian = read_w90_tb(dataset"Si2_valence/outputs/WS/Si2_valence").hamiltonian
    hamiltonian = read_w90_hr(
        dataset"Si2_valence/outputs/WS/Si2_valence", ref_hamiltonian.Rspace.lattice
    )
    @test isapprox(hamiltonian, ref_hamiltonian; atol = 1.0e-5)
end

@testitem "read_w90_hr MDRS" begin
    using Wannier.Datasets
    ref_hamiltonian = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence").hamiltonian
    hamiltonian = read_w90_hr(
        dataset"Si2_valence/outputs/MDRS/Si2_valence", ref_hamiltonian.Rspace.lattice
    )
    @test isapprox(hamiltonian, ref_hamiltonian; atol = 1.0e-5)
end

@testitem "write_w90_hr" begin
    using Wannier.Datasets
    ref_hamiltonian = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence").hamiltonian

    outdir = mktempdir(; cleanup = true)
    outprefix = joinpath(outdir, "Si2_valence")
    write_w90_hr(outprefix, ref_hamiltonian)

    hamiltonian = read_w90_hr(outprefix, ref_hamiltonian.Rspace.lattice)
    @test isapprox(ref_hamiltonian, hamiltonian; atol = 1.0e-5)
end

@testitem "read_w90_tb_chk_spn" begin
    using Wannier.Datasets

    hamiltonian, position, spin = Wannier.read_w90_tb_chk_spn(
        dataset"Fe_soc_coarse/outputs/Fe";
        chk = dataset"Fe_soc_coarse/outputs/Fe.chk",
        spn = dataset"Fe_soc_coarse/Fe.spn",
    )

    S11 = [
        0.03333792819623771 - 2.9369081252907633e-18im
        -0.04529527445929681 - 1.5064650151202397e-18im
        0.046510020098334805 + 1.7122606109611054e-18im
    ]
    S12 = [
        -0.03221828528841446 + 0.00615582634310275im
        0.0014468584576766712 - 0.021998294535299574im
        -0.02024967831556668 - 0.04452681873304288im
    ]
    S21 = [0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]
    S22 = [
        0.032613311617921685 - 1.0347471943305515e-17im
        0.008923620770530381 - 1.1770746615814775e-17im
        0.04248374526138618 + 1.3246640315040066e-17im
    ]
    ref_S = [[S11, S21] [S12, S22]]
    @test spin[1][1:2, 1:2] ≈ ref_S
end
