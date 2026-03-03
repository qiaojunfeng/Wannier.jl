module WannierMakieExt

using Makie

using Wannier
import Wannier: bandplot, bandplot!, get_bandplot

@recipe BandPlot (kpath, eigenvals) begin
    "Fermi energy"
    fermi_energy = nothing

    "Whether to shift the bands such that the Fermi energy is at zero"
    shift_fermi = false

    "Legend label for the bands"
    label = nothing

    "Line color"
    linecolor = @inherit linecolor :black

    "Line style"
    linestyle = @inherit linestyle :solid

    "Line width"
    linewidth = @inherit linewidth 1

    "Line color for vertical lines at high-symmetry kpoints"
    linecolor_ksym = :grey

    "Line style for vertical lines"
    linestyle_ksym = :solid

    "Vertical line width"
    linewidth_ksym = 0.5

    "Whether to show horizontal line at Fermi energy"
    show_fermi = true

    "Line color for horizontal line at Fermi energy"
    linecolor_fermi = :blue

    "Line style for horizontal line at Fermi energy"
    linestyle_fermi = :dash

    "Line width for horizontal line at Fermi energy"
    linewidth_fermi = 0.8

    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(
    p::BandPlot{Tuple{K,E}}
) where {K<:Wannier.RecipPath,E<:AbstractVector{<:AbstractVector{<:Real}}}
    # Cumulative distance along kpath
    x = Wannier.get_linear_path(p.kpath[])
    # Eigenvalues of the bands
    eigenvals = p.eigenvals[]

    nkpts = length(eigenvals)
    nkpts == length(x) || error("Length of eigenvals does not match length of kpoints")
    nbands = length(eigenvals[1])

    # Labels of high-symmetry kpoints
    labs = Wannier.symm_label_to_unicode(p.kpath[].labels)
    # Indices and labels of merged high-symmetry kpoints
    xt_idxs, xt_labs = Wannier.merge_symm_labels(p.kpath[].indices, labs)
    # Actual x positions of the high-symmetry kpoints
    xts = x[xt_idxs]

    vlines!(
        p,
        xts;
        color=p.linecolor_ksym[],
        linewidth=p.linewidth_ksym[],
        linestyle=p.linestyle_ksym[],
    )

    map!(p.attributes, [:fermi_energy, :shift_fermi], :_fermi_energy) do fermi_energy, shift_fermi
        if isnothing(fermi_energy)
            return nothing
        else
            return shift_fermi ? 0.0 : fermi_energy
        end
    end
    if p.show_fermi[] && !isnothing(p._fermi_energy[])
        hlines!(
            p,
            p._fermi_energy[];
            color=p.linecolor_fermi[],
            linewidth=p.linewidth_fermi[],
            linestyle=p.linestyle_fermi[],
        )
    end

    map!(p.attributes, [:fermi_energy, :shift_fermi], :bands) do fermi_energy, shift_fermi
        # The input eigenvals is a vector of vectors, where the outer vector is
        # over k-points and the inner vector is over bands. We need to reshape
        # it to a vector of length nbands with inner vector of length nkpoints
        # for plotting.
        bands = map(1:nbands) do ib
            [eigenvals[ik][ib] for ik in 1:nkpts]
        end
        if shift_fermi && !isnothing(fermi_energy)
            return map(x -> x .- fermi_energy, bands)
        else
            return bands
        end
    end

    for (i, bs) in enumerate(p.bands[])
        kwargs = (; color=p.linecolor[], linestyle=p.linestyle[], linewidth=p.linewidth[])
        # Only show label for the first band, to avoid cluttering the legend
        (i == 1) && (kwargs = (; kwargs..., label=p.label[]))
        lines!(p, x, bs; kwargs...)
    end

    return p
end

function fig_ax_bandplot(kpath::Wannier.RecipPath; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])

    shift_fermi = get(kwargs, :shift_fermi, false)
    ylabel = get(kwargs, :ylabel, nothing)
    if isnothing(ylabel)
        if shift_fermi
            ylabel = rich(rich("E - E"; font=:italic), subscript("F"), " (eV)")
        else
            ylabel = "Energy (eV)"
        end
    end
    ax.ylabel = ylabel

    # Cumulative distance along kpath
    x = Wannier.get_linear_path(kpath)
    # Labels of high-symmetry kpoints
    labs = Wannier.symm_label_to_unicode(kpath.labels)
    # Indices and labels of merged high-symmetry kpoints
    xt_idxs, xt_labs = Wannier.merge_symm_labels(kpath.indices, labs)
    # Actual x positions of the high-symmetry kpoints
    xts = x[xt_idxs]
    ax.xticks = (xts, xt_labs)

    ax.ygridvisible = false

    xlims!(ax, [minimum(x), maximum(x)])

    return fig, ax
end

#=
This is a bit hacky: the recipe only works with the Plot, but we want to
set the ylabel of the Axis, etc. Therefore, we need this wrapper function.
=#
function Wannier.get_bandplot(
    kpath::Wannier.RecipPath, eigenvals::E; kwargs...
) where {E<:AbstractVector{<:AbstractVector{<:Real}}}
    fig, ax = fig_ax_bandplot(kpath; kwargs...)
    p = bandplot!(ax, kpath, eigenvals; kwargs...)
    return Makie.FigureAxisPlot(fig, ax, p)
end

function Wannier.get_bandplot(
    kpath::Wannier.RecipPath, eigenvals1::E, eigenvals2::E; kwargs1=(;), kwargs2=(;)
) where {E<:AbstractVector{<:AbstractVector{<:Real}}}
    fig, ax = fig_ax_bandplot(kpath; kwargs1...)

    color1 = get(kwargs1, :linecolor, :black)
    color2 = get(kwargs2, :linecolor, :red)
    #
    linestyle1 = get(kwargs1, :linestyle, :solid)
    linestyle2 = get(kwargs2, :linestyle, :dash)
    #
    label1 = get(kwargs1, :label, "Band 1")
    label2 = get(kwargs2, :label, "Band 2")
    #
    kwargs1 = (; kwargs1..., linecolor=color1, linestyle=linestyle1, label=label1)
    kwargs2 = (; kwargs2..., linecolor=color2, linestyle=linestyle2, label=label2)

    plots = [
        bandplot!(ax, kpath, eigenvals1; kwargs1...)
        bandplot!(ax, kpath, eigenvals2; kwargs2...)
    ]

    axislegend(ax)

    # return (; fig, ax, plots)
    # Just return the primary plot, as with FigureAxisPlot, the returned plot
    # is automatically displaed.
    return Makie.FigureAxisPlot(fig, ax, plots[1])
end

end # module
