@testitem "BerryCurvatureInterpolator" begin
    using LinearAlgebra
    using DelimitedFiles
    using Wannier.Datasets
    hamiltonian, position = read_w90_tb(dataset"Fe/reference/MDRS/Fe")
    interp = Wannier.BerryCurvatureInterpolator(hamiltonian, position)

    ref_kpt = read_w90_band_kpt(dataset"Fe/reference/MDRS/postw90/Fe-path.kpt")
    ref_dat = readdlm(dataset"Fe/reference/MDRS/postw90/Fe-curv.dat")
    # w90 actually writes -Ω, so we need to negate it
    ref_Ω = map(eachrow(ref_dat[:, 2:end])) do Ω
        -Wannier.vector_to_tensor(Ω)
    end

    # if I use the kpoints in ref_kpt, the difference between eigenvalues is
    # around 1e-4, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_kpt.kpoints
    win = read_win(dataset"Fe/Fe.win")
    kpi = generate_w90_kpoint_path(win.unit_cell_cart, win.kpoint_path)
    # postw90.x has a bug, it misses the `H` point at 417
    kpoints = get_kpoints(kpi)
    deleteat!(kpoints, 417)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1e-6)

    interp_H = HamiltonianInterpolator(hamiltonian)
    eigenvalues = interp_H(kpoints)[1]
    occupations = [Int.(e .<= win.fermi_energy) for e in eigenvalues]
    Ω = interp(kpoints, occupations)

    Ω_band = interp(kpoints)
    Ω_band_sum = map(zip(Ω_band, occupations)) do (Ωₖ, fₖ)
        sum(fₖ .* Ωₖ)
    end
    # the Ω is inevitably more accurate than Ω_band_sum, so we use a rather
    # loose tolerance here
    @assert all(isapprox.(Ω, Ω_band_sum; atol=5e-2))

    # TODO understand why they are different
    @test all(isapprox.(Ω, ref_Ω; atol=2e-1))
end
