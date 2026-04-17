@testmodule OrbitalMagnetizationEnv begin
    using Wannier
    using Wannier.Datasets
    export model, hamiltonian_position, position_hamiltonian_position, uHu

    model = read_w90_with_chk(dataset"Fe_soc_coarse/Fe", dataset"Fe_soc_coarse/outputs/Fe.chk")
    hamiltonian_position = TBHamiltonianPosition(model)
    uHu = read_uHu(dataset"Fe_soc_coarse/Fe.uHu").uHu
    position_hamiltonian_position = TBPositionHamiltonianPosition(model, uHu)
end

@testitem "TBHamiltonianPosition" setup = [OrbitalMagnetizationEnv] begin
    using LinearAlgebra

    Hr_11 = Vec3(
        -0.1310914834482943 + 0.05853108710468198im,
        -0.10609815871369616 - 0.06165462741261124im,
        -0.22750100392145212 + 0.12742014536973748im
    )
    Hr_12 = Vec3(
        0.0 + 0.0im,
        0.0 + 0.0im,
        0.0 + 0.0im
    )
    Hr_21 = Vec3(
        0.2858787836391398 + 0.5800944404791165im,
        0.07700692054319432 - 0.24674019231253563im,
        -0.6998327066540999 - 0.6544517149242701im
    )
    Hr_22 = Vec3(
        -0.08208690744174643 + 0.122948638299064im,
        -0.13302713637817998 - 0.13666634789988658im,
        -0.004434708333426407 - 0.12706594346508193im
    )
    ref_Hr_001 = [[Hr_11, Hr_21] [Hr_12, Hr_22]]
    @test isapprox(hamiltonian_position[0, 1, 0][1:2, 1:2], ref_Hr_001; atol = 1.0e-8)
end

@testitem "TBPositionHamiltonianPosition" setup = [OrbitalMagnetizationEnv] begin
    using LinearAlgebra

    ref_rHr_011 = [
        [
            -0.011344631721247528 - 8.797323657074556e-18im   0.08061997095049099 + 0.055782599683839265im       0.1512117475770554 - 0.008811001411659518im
            0.08061997095049099 - 0.055782599683839265im   -0.06864220679169453 + 1.3972282670756823e-17im    -0.1294603039919712 + 0.08477139042404999im
            0.1512117475770554 + 0.008811001411659518im    -0.1294603039919712 - 0.08477139042404999im     -0.024328266922674535 + 1.6971365538098622e-17im
        ],
        [
            -0.13102251317841107 - 0.0872578912067857im     0.22701200284340356 + 0.0552532493128975im    0.27189683515426977 + 0.21875009920028324im
            0.3298083924817641 + 0.16327194339453072im     -0.076738703576335 - 0.010582808772552924im  -0.3573711109596745 - 0.24639555412040082im
            0.24099374405778362 + 0.056959474467350446im  -0.21820206700118655 + 0.03269998910832419im   0.10510780019703406 + 0.08897746731961312im
        ],
    ]
    @test isapprox(position_hamiltonian_position[0, 1, 0][1:2, 1], ref_rHr_011; atol = 1.0e-8)
end

@testitem "OrbitalMagnetizationInterpolator" setup = [OrbitalMagnetizationEnv] begin
    using LinearAlgebra
    using DelimitedFiles
    using Wannier.Datasets

    # note that when w90 writes tb.dat, it use imag(log(...)) for
    # the diagonal part of position operator, thus it will be different from
    # the one used in postw90.x, which directly use the overlap matrices for
    # position operator. Therefore, we read directly from chk file to reproduce
    # the same results.
    # hamiltonian, position = read_w90_tb(dataset"Fe_soc_coarse/outputs/Fe")
    hamiltonian = TBHamiltonian(model)
    Rspace = generate_Rspace(model)
    position = TBPosition(Rspace, model; imlog_diag = false)
    hamiltonian_position1 = TBHamiltonianPosition(Rspace, model)
    position_hamiltonian_position1 = TBPositionHamiltonianPosition(Rspace, model, uHu)
    win = read_win(dataset"Fe_soc_coarse/Fe.win")
    interp = Wannier.OrbitalMagnetizationInterpolator(
        hamiltonian,
        position,
        hamiltonian_position1,
        position_hamiltonian_position1,
        win["fermi_energy"],
    )

    ref_kpt = read_w90_band_kpt(dataset"Fe_soc_coarse/outputs/postw90/Fe-path.kpt")
    ref_dat = readdlm(dataset"Fe_soc_coarse/outputs/postw90/Fe-morb.dat")
    # w90 actually writes -1/2 * M, where M = LVTS12 Eq. 97
    ref_M = map(eachrow(ref_dat[:, 2:end])) do M
        -2 * Wannier.axialvector_to_antisymmetrictensor(M)
    end

    # if I use the kpoints in ref_kpt, the difference between eigenvalues is
    # around 1e-4, this is because the kpoints coordinates do not have enough
    # digits. Therefore, I read the win file and construct the kpoints myself.
    # kpoints = ref_kpt.kpoints
    kseg = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
    # I use less kpoints to speedup the test
    kpath = KPath(kseg, 5)
    kpoints = collect(kpath)
    # postw90.x has a bug, it misses the `H` point
    ik_H = 22
    deleteat!(kpoints, ik_H)
    @test all(norm.(kpoints - ref_kpt.kpoints) .< 1.0e-6)

    M = interp(kpoints)
    @test all(isapprox.(M, ref_M; atol = 5.0e-7))
end
