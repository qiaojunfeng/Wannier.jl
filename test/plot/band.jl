@testitem "bandplot" begin
    using Wannier.Datasets
    using CairoMakie

    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = reciprocal_lattice(win["unit_cell_cart"])
    kpath, eigenvals = read_w90_band(dataset"Si2_valence/outputs/MDRS/Si2_valence", recip_lattice)

    # fig, ax, p = bandplot(kpath, eigenvals; label="Wan")

    δ = 0.2
    eigenvals2 = collect(eachcol(reduce(hcat, eigenvals) .+ δ))
    kwargs1 = (; label = "Wan1")
    kwargs2 = (; label = "Wan2")
    fig, ax, p = Wannier.get_bandplot(
        kpath,
        eigenvals,
        eigenvals2;
        kwargs1,
        kwargs2,
        fermi_energy = win["fermi_energy"],
        shift_fermi = true,
    )

    @test all(x -> x isa Plot{Wannier.bandplot}, ax.scene.plots)
    @test p.attributes[:fermi_energy][] == win["fermi_energy"]
end

@testitem "projbandplot" begin
    using Wannier.Datasets
    using CairoMakie

    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = reciprocal_lattice(win["unit_cell_cart"])
    kpath, eigenvals = read_w90_band(dataset"Si2_valence/outputs/MDRS/Si2_valence", recip_lattice)
    # Fake projections
    nprojs = 2
    U = rand_gauge(ComplexF64, length(kpath), length(eigenvals[1]), nprojs)
    projs = Wannier.projectability(U)
    labels = ["WF $i" for i in 1:nprojs]

    kwargs = (; fermi_energy = win["fermi_energy"], shift_fermi = true)
    fig, ax, p = Wannier.get_projbandplot(kpath, eigenvals, projs, labels; kwargs...)

    @test all(x -> x isa Plot{Wannier.projbandplot}, ax.scene.plots)
    @test p.attributes[:fermi_energy][] == win["fermi_energy"]
end
