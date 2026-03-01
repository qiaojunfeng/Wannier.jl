@testitem "Hamiltonian MDRS" begin
    using Wannier.Datasets
    using CairoMakie

    win = read_win(dataset"Si2_valence/Si2_valence.win")
    recip_lattice = reciprocal_lattice(win.unit_cell_cart)
    kpi, eigenvals = read_w90_band(dataset"Si2_valence/outputs/MDRS/Si2_valence", recip_lattice)
    kpath = Wannier.RecipPath(kpi)

    fig, ax, p = bandplot(kpath, eigenvals; label="Wan")

    @test any(x -> x isa BandPlot, ax.scene.plots)

    @test p.attributes[:color][] == :red
end
