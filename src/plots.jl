import Plots as Pl
using .Parameters: Bands, Projectabilities

function plot_bands(bands::Bands)
    Pl.PlotlyBackend()
    plt = Pl.current()
    return plot_bands(plt, bands)
end

function plot_bands(plt, bands::Bands; fermi_energy::Union{Int,Float64}=0.0, 
    colors::Union{Matrix{Float64},Missing}=missing)
# t = range(0, stop = 1, length = 100)
# θ = (6π) .* t
# x = t .* cos.(θ)
# y = t .* sin.(θ)
# p1 = plot(x, y, line_z = t, linewidth = 3, legend = false)
# p2 = scatter(x, y, marker_z = (+), color = :bluesreds, legend = false)
# plot(p1, p2)
    Pl.ylabel!(plt, "Energy (eV)")

    Pl.hline!(plt, [fermi_energy]; linestyle=:dash, linecolor=:blue, linewidth=0.2)

    if bands.num_symm_points != 0
        Pl.vline!(plt, bands.kpaths[bands.symm_points]; linecolor=:black, linewidth=0.2)
        Pl.xticks!(plt, bands.kpaths[bands.symm_points], bands.symm_points_label)
    end

    if isa(colors, Missing)
        Pl.plot!(plt, bands.kpaths, bands.energies; legend=false, grid=false, framestyle=:box)
    else
        Pl.plot!(plt, bands.kpaths, bands.energies, line_z=colors; 
        color=:copper, legend=false, grid=false, framestyle=:box)
    end

    Pl.xlims!(plt, (bands.kpaths[1], bands.kpaths[end]))

    emin = minimum(bands.energies) - 0.5
    emax = maximum(bands.energies) + 0.5
    Pl.ylims!(plt, (emin, emax))

    return plt
end

function plot_bands_projectabilities(bands::Bands, projectabilities::Projectabilities; 
    fermi_energy::Union{Int,Float64}=0.0)

    colors = dropdims(sum(projectabilities.proj, dims=3), dims=3)

    plt = plot_bands(Pl.plot(), bands; fermi_energy=fermi_energy, colors=colors)

    Pl.gui(plt)

    return plt
end

function plot_bands_diff(bands1::Bands, bands2::Bands; fermi_energy::Union{Int,Float64}=0.0, 
    label1::String="Bands1", label2::String="Bands2")
    plt = Pl.plot(grid=false, framestyle=:box)

    Pl.ylabel!(plt, "Energy (eV)")
    Pl.hline!(plt, [fermi_energy]; linestyle=:dash, linecolor=:blue, linewidth=0.2, label="Fermi energy")

    if bands1.num_symm_points != 0 && all(bands1.symm_points_label .!= "")
        Pl.vline!(plt, bands1.kpaths[bands1.symm_points]; linecolor=:black, linewidth=0.2, label="")
        Pl.xticks!(plt, bands1.kpaths[bands1.symm_points], bands1.symm_points_label)
    end
    if bands2.num_symm_points != 0 && all(bands2.symm_points_label .!= "")
        Pl.vline!(plt, bands2.kpaths[bands2.symm_points]; linecolor=:black, linewidth=0.2, label="")
        Pl.xticks!(plt, bands2.kpaths[bands2.symm_points], bands2.symm_points_label)
    end

    # only show 1 label in the legend
    labels = hcat(label1, fill("", 1, bands1.num_bands))
    Pl.plot!(plt, bands1.kpaths, bands1.energies; linecolor=:black, label=labels)
    labels = hcat(label2, fill("", 1, bands2.num_bands))
    Pl.plot!(plt, bands2.kpaths, bands2.energies; linecolor=:red, linestyle=:dash, label=labels)

    xmin = min(bands1.kpaths[1], bands2.kpaths[1])
    xmax = max(bands1.kpaths[end], bands2.kpaths[end])
    Pl.xlims!(plt, (xmin, xmax))

    emin = min(minimum(bands1.energies), minimum(bands2.energies)) - 0.5
    emax = max(maximum(bands1.energies), maximum(bands2.energies)) + 0.5
    Pl.ylims!(plt, (emin, emax))

    Pl.gui(plt)

    return plt
end
