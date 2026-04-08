@testitem "BerryCurvatureInterpolator" begin
    using LinearAlgebra
    using DelimitedFiles
    using Wannier.Datasets
    # note that when w90 writes tb.dat, it use imag(log(...)) for
    # the diagonal part of position operator, thus it will be different from
    # the one used in postw90.x, which directly use the overlap matrices for
    # position operator. Therefore, we read directly from chk file to reproduce
    # the same results.
    # hamiltonian, position = read_w90_tb(dataset"Fe_soc/outputs/MDRS/Fe")
    model = read_w90_with_chk(dataset"Fe_soc/Fe", dataset"Fe_soc/outputs/Fe.chk")
    hamiltonian = TBHamiltonian(model)
    position = TBPosition(model; imlog_diag = false)
    win = read_win(dataset"Fe_soc/Fe.win")
    interp = Wannier.BerryCurvatureInterpolator(hamiltonian, position, win["fermi_energy"])

    # I use less kpoints to speedup the test
    COARSE = true
    postw90_dir = COARSE ? "postw90_coarse" : "postw90"

    ref_kpt = read_w90_band_kpt(dataset"Fe_soc/outputs/MDRS/" * postw90_dir * "/Fe-path.kpt")
    ref_dat = readdlm(dataset"Fe_soc/outputs/MDRS/" * postw90_dir * "/Fe-curv.dat")
    # w90 actually writes -Ω, so we need to negate it
    ref_Ω = map(eachrow(ref_dat[:, 2:end])) do Ω
        -Wannier.axialvector_to_antisymmetrictensor(Ω)
    end

    # if I use the kpoints in ref_kpt, the difference between eigenvalues is
    # around 1e-4, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_kpt.kpoints
    kseg = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    kpath = KPath(kseg, COARSE ? 5 : 100)
    kpoints = collect(kpath)
    # postw90.x has a bug, it misses the `H` point
    ik_H = COARSE ? 22 : 417
    deleteat!(kpoints, ik_H)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1.0e-6)

    # summed over bands
    Ω = interp(kpoints, Wannier.WYSV06())
    @test all(isapprox.(Ω, ref_Ω; atol = 5.0e-6))

    # band-resolved Berry curvature
    Ω_band = interp(kpoints, Wannier.WYSV06BandResolved())
    eigenvalues = HamiltonianInterpolator(hamiltonian)(kpoints)[1]
    occupations = [Int.(εₖ .<= win["fermi_energy"]) for εₖ in eigenvalues]
    Ω_band_sum = map(zip(Ω_band, occupations)) do (Ωₖ, fₖ)
        sum(fₖ .* Ωₖ)
    end
    @test all(isapprox.(Ω, Ω_band_sum; atol = 1.0e-10))

    Ω = interp(kpoints, Wannier.LVTS12())
    @test all(isapprox.(Ω, ref_Ω; atol = 5.0e-6))
end
