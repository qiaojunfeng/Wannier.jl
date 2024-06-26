@testitem "spin projection" begin
    using LinearAlgebra
    using Wannier.Datasets
    hamiltonian, position, spin = read_w90_tb_chk_spn(
        dataset"Fe_soc/reference/MDRS/Fe";
        spn=dataset"Fe_soc/Fe.spn",
        chk=dataset"Fe_soc/reference/Fe.chk",
    )
    # project onto the z axis
    θ = 0.0
    ϕ = 0.0
    interp = SpinProjectionInterpolator(hamiltonian, spin, θ, ϕ)

    ref_kpt = read_w90_band_kpt(dataset"Fe_soc/reference/MDRS/postw90/Fe-path.kpt")
    ref_dat = read_w90_band_dat(dataset"Fe_soc/reference/MDRS/postw90/Fe-bands.dat")

    # if I use the kpoints in ref_kpt, the difference between eigenvalues is
    # around 1e-4, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_kpt.kpoints
    win = read_win(dataset"Fe_soc/Fe.win")
    kpi = generate_w90_kpoint_path(win.unit_cell_cart, win.kpoint_path)
    # postw90.x has a bug, it misses the `H` point at 417
    kpoints = get_kpoints(kpi)
    deleteat!(kpoints, 417)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1e-6)
    ##
    eigenvalues = HamiltonianInterpolator(hamiltonian)(kpoints)[1]
    @test all(norm.(eigenvalues - ref_dat.eigenvalues) .< 2e-6)

    Sz = interp(kpoints)
    @test all(isapprox.(Sz, ref_dat.extras; atol=5e-5))
end
