@testitem "bandplot" begin
    using Wannier.Datasets
    using CairoMakie

    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = reciprocal_lattice(win.unit_cell_cart)
    kpi, eigenvals = read_w90_band(dataset"Si2_valence/outputs/MDRS/Si2_valence", recip_lattice)
    kpath = Wannier.RecipPath(kpi)

    # fig, ax, p = bandplot(kpath, eigenvals; label="Wan")

    δ = 0.2
    eigenvals2 = map(x->x .+ δ, eigenvals)
    kwargs1 = (; label="Wan1", win.fermi_energy, shift_fermi=true)
    kwargs2 = (; kwargs1..., label="Wan2")
    fig, ax, p = Wannier.get_bandplot(kpath, eigenvals, eigenvals2; kwargs1, kwargs2)

    @test all(x -> x isa Plot{Wannier.bandplot}, ax.scene.plots)
    @test p.attributes[:fermi_energy][] == win.fermi_energy
end

@testitem "projbandplot" begin
    using Wannier.Datasets
    using CairoMakie

    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = reciprocal_lattice(win.unit_cell_cart)
    kpi, eigenvals = read_w90_band(dataset"Si2_valence/outputs/MDRS/Si2_valence", recip_lattice)
    kpath = Wannier.RecipPath(kpi)
    # Fake projections
    nprojs = 2
    projections = rand_gauge(ComplexF64, length(kpi), length(eigenvals[1]), nprojs)
    projections = Wannier.compute_projectability(projections, [[i] for i in 1:nprojs])
    labels = ["WF $i" for i in 1:nprojs]

    kwargs = (; win.fermi_energy, shift_fermi=true)
    fig, ax, p = Wannier.get_projbandplot(kpath, eigenvals, projections, labels; kwargs...)

    @test all(x -> x isa Plot{Wannier.projbandplot}, ax.scene.plots)
    @test p.attributes[:fermi_energy][] == win.fermi_energy
end
