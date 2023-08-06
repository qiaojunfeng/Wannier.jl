@testitem "TBHamiltonianPosition" begin
    using LinearAlgebra
    using Wannier: Vec3
    using Wannier.Datasets

    model = read_w90_with_chk(dataset"Fe_soc/Fe", dataset"Fe_soc/reference/Fe.chk")
    hamiltonian_position = TBHamiltonianPosition(model)

    Hr_11 = Vec3(
        -0.05655034287285573 - 0.0033676997474095015im,
        0.14573076408128938 + 0.002481532972539341im,
        -0.08495167782844368 + 0.0002025296058253332im,
    )
    Hr_12 = Vec3(
        -0.004894341665464354 - 0.0035083746559635472im,
        8.947905061703487e-5 - 0.0003222480961347094im,
        -0.002246877325241364 - 0.010270224128882639im,
    )
    Hr_21 = Vec3(
        -0.008879207283281055 + 0.017350636753024708im,
        -0.001550320638839293 + 0.011450968898925317im,
        -0.008105982416168542 + 0.0016898640836712176im,
    )
    Hr_22 = Vec3(
        -0.0420601683301116 - 0.00033243732755332066im,
        0.08705831014849093 + 0.009332448985531172im,
        0.17436896103376767 - 0.0007776111447994936im,
    )
    ref_Hr_001 = [[Hr_11, Hr_21] [Hr_12, Hr_22]]
    @test isapprox(hamiltonian_position[0, 1, 1][1:2, 1:2], ref_Hr_001; atol=1e-8)
end

@testitem "TBPositionHamiltonianPosition" begin
    using LinearAlgebra
    using Wannier: Vec3
    using Wannier.Datasets

    model = read_w90_with_chk(dataset"Fe_soc/Fe", dataset"Fe_soc/reference/Fe.chk")
    uHu = read_uHu(dataset"Fe_soc/Fe.uHu")
    position_hamiltonian_position = TBPositionHamiltonianPosition(model, uHu)

    ref_rHr_011 = [
        [
            0.1791744538931292+0.0020955816183890405im -0.186068778524343-0.0020362406364252553im -0.09391057298063965-0.0018317497603643198im
            0.04124275653175067-0.009008304363442845im 0.04282298381236467-0.0029005984505842833im -0.13165388156091815+0.008739659094908092im
            -0.24195275404031225-0.005679662165088098im 0.13953004528724272+0.0035418448818636747im 0.20921953895157094+0.006034632637907919im
        ],
        [
            0.03957686145822008-0.06685710036635506im -0.003926564560497527+0.0023361688298376535im 0.008057284311593447-0.0030866754673423806im
            0.001472662662175219+0.0006408293375997858im -0.014179847648214864+0.03323642945856213im -0.0036383288813261407-0.002937128744984542im
            -0.0027410446394830847+0.0068988661552413285im 0.002527434016376507+0.0033246397145206582im -0.002834534631742737-0.0072179750352187625im
        ],
    ]
    @test isapprox(position_hamiltonian_position[0, 1, 1][1:2, 1], ref_rHr_011; atol=1e-8)
end

@testitem "OrbitalMagnetizationInterpolator" begin
    using LinearAlgebra
    using DelimitedFiles
    using Wannier.Datasets
    # note that when w90 writes tb.dat, it use imag(log(...)) for
    # the diagonal part of position operator, thus it will be different from
    # the one used in postw90.x, which directly use the overlap matrices for
    # position operator. Therefore, we read directly from chk file to reproduce
    # the same results.
    # hamiltonian, position = read_w90_tb(dataset"Fe_soc/reference/MDRS/Fe")
    model = read_w90_with_chk(dataset"Fe_soc/Fe", dataset"Fe_soc/reference/Fe.chk")
    hamiltonian = TBHamiltonian(model)
    Rspace = generate_Rspace(model)
    position = TBPosition(Rspace, model; imlog_diag=false)
    hamiltonian_position = TBHamiltonianPosition(Rspace, model)
    uHu = read_uHu(dataset"Fe_soc/Fe.uHu")
    position_hamiltonian_position = TBPositionHamiltonianPosition(Rspace, model, uHu)
    win = read_win(dataset"Fe_soc/Fe.win")
    interp = Wannier.OrbitalMagnetizationInterpolator(
        hamiltonian,
        position,
        hamiltonian_position,
        position_hamiltonian_position,
        win.fermi_energy,
    )

    ref_kpt = read_w90_band_kpt(dataset"Fe_soc/reference/MDRS/postw90/Fe-path.kpt")
    ref_dat = readdlm(dataset"Fe_soc/reference/MDRS/postw90/Fe-morb.dat")
    # w90 actually writes -1/2 * M, where M = LVTS12 Eq. 97
    ref_M = map(eachrow(ref_dat[:, 2:end])) do M
        -2 * Wannier.axialvector_to_antisymmetrictensor(M)
    end

    # if I use the kpoints in ref_kpt, the difference between eigenvalues is
    # around 1e-4, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_kpt.kpoints
    kpi = generate_w90_kpoint_path(win.unit_cell_cart, win.kpoint_path)
    # postw90.x has a bug, it misses the `H` point at 417
    kpoints = get_kpoints(kpi)
    deleteat!(kpoints, 417)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1e-6)

    M = interp(kpoints)
    @test all(isapprox.(M, ref_M; atol=5e-7))
end
