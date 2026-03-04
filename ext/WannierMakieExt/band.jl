import Wannier: bandplot, bandplot!, get_bandplot

function shared_bandplot_attributes()
    Makie.@DocumentedAttributes begin
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
    end
end

@recipe BandPlot (kpath, eigenvals) begin
    shared_bandplot_attributes()...

    Makie.mixin_generic_plot_attributes()...
end

function get_xtick_indices_labels(kpath::Wannier.RecipPath)
    # Labels of high-symmetry kpoints
    labs = Wannier.symm_label_to_unicode(kpath.labels)
    # Indices and labels of merged high-symmetry kpoints
    xt_idxs, xt_labs = Wannier.merge_symm_labels(kpath.indices, labs)
    return xt_idxs, xt_labs
end

function Makie.plot!(
    p::BandPlot{Tuple{K,E}}
) where {K<:Wannier.RecipPath,E<:AbstractVector{<:AbstractVector{<:Real}}}
    map!(p.attributes, [:kpath], :nkpts) do kpath
        return length(kpath)
    end
    map!(p.attributes, [:eigenvals], :nbands) do eigenvals
        p.nkpts[] == length(eigenvals) || error("Length of eigenvals does not match length of kpoints")
        return length(eigenvals[1])
    end
    # Cumulative distance along kpath
    map!(p.attributes, [:kpath], :x) do kpath
        return Wannier.get_linear_path(kpath)
    end
    # Positions of vertical lines at high-symmetry kpoints
    map!(p.attributes, [:kpath, :x], :x_ksym) do kpath, x
        xt_idxs = get_xtick_indices_labels(kpath)[1]
        return x[xt_idxs]
    end
    if !isempty(p.x_ksym[])
        vlines!(
            p,
            p.x_ksym[];
            color=p.linecolor_ksym[],
            linewidth=p.linewidth_ksym[],
            linestyle=p.linestyle_ksym[],
        )
    end

    map!(
        p.attributes, [:fermi_energy, :shift_fermi], :y_fermi
    ) do fermi_energy, shift_fermi
        if isnothing(fermi_energy)
            return nothing
        else
            return shift_fermi ? 0.0 : fermi_energy
        end
    end
    if p.show_fermi[] && !isnothing(p.y_fermi[])
        hlines!(
            p,
            p.y_fermi[];
            color=p.linecolor_fermi[],
            linewidth=p.linewidth_fermi[],
            linestyle=p.linestyle_fermi[],
        )
    end

    map!(p.attributes, [:fermi_energy, :shift_fermi], :bands) do fermi_energy, shift_fermi
        # The input eigenvals is a vector of vectors, where the outer vector is
        # over k-points and the inner vector is over bands. We need to reshape
        # it such that it is accessed by bands[ib][ik].
        bands = map(1:p.nbands[]) do ib
            [p.eigenvals[][ik][ib] for ik in 1:p.nkpts[]]
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
        lines!(p, p.x[], bs; kwargs...)
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
            # Note below I am using the unicode minus sign (U+2212) instead of
            # the ASCII hyphen (U+002D)
            ylabel = rich(rich("E − E"; font=:italic), subscript("F"), " (eV)")
        else
            ylabel = "Energy (eV)"
        end
    end
    ax.ylabel = ylabel
    ax.ygridvisible = false

    x = Wannier.get_linear_path(kpath)
    xlims!(ax, [minimum(x), maximum(x)])

    xt_idxs, xt_labs = get_xtick_indices_labels(kpath)
    # Actual x positions of the high-symmetry kpoints
    xt = x[xt_idxs]
    ax.xticks = (xt, xt_labs)

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
