@testitem "generate_Rspace_domain WS" begin
    using Wannier.Datasets
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    tbdat = read_w90_tbdat(dataset"Si2_valence/reference/ws/Si2_valence_tb.dat")

    ref_Rdomain = Wannier.WSRspaceDomain(tbdat.lattice, tbdat.Rvectors, tbdat.Rdegens)
    Rdomain = generate_Rspace_domain(
        ref_Rdomain.lattice, win.mp_grid, Wannier.WSInterpolation()
    )
    @test Rdomain ≈ ref_Rdomain
end

@testitem "generate_Rspace_domain MDRS" begin
    using Wannier.Datasets
    win = read_win(dataset"Si2_valence/Si2_valence.win")
    wout = read_wout(dataset"Si2_valence/reference/Si2_valence.wout")
    wsvec = read_w90_wsvec(dataset"Si2_valence/reference/mdrs/Si2_valence_wsvec.dat")
    tbdat = read_w90_tbdat(dataset"Si2_valence/reference/mdrs/Si2_valence_tb.dat")

    ref_Rdomain = Wannier.MDRSRspaceDomain(
        tbdat.lattice, tbdat.Rvectors, tbdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens
    )
    # to fractional
    inv_lattice = inv(tbdat.lattice)
    centers = map(wout.centers) do c
        inv_lattice * c
    end
    Rdomain = generate_Rspace_domain(
        ref_Rdomain.lattice, win.mp_grid, centers, Wannier.MDRSInterpolation()
    )
    @test Rdomain ≈ ref_Rdomain
end
