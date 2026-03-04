import Wannier: projbandplot, projbandplot!, get_projbandplot

@recipe ProjBandPlot (kpath, eigenvals, projs, labels) begin
    shared_bandplot_attributes()...

    "Marker size for the scatter points representing the projectabilities"
    markersize = @inherit markersize 10

    "Colors for the orbitals. If not specified, will be taken from the color scheme."
    colors = nothing

    "Color scheme for the orbitals"
    colorscheme = :Set1_9
    # colorscheme = :tab10

    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(
    p::ProjBandPlot{Tuple{K,E,P,L}}
) where {
    K<:Wannier.RecipPath,
    E<:AbstractVector{<:AbstractVector{<:Real}},
    P<:AbstractVector{<:AbstractMatrix{<:Real}},
    L<:AbstractVector{<:AbstractString},
}
    map!(p.attributes, [:eigenvals], :nkpts) do eigenvals
        return length(eigenvals)
    end
    map!(p.attributes, [:eigenvals], :nbands) do eigenvals
        return length(eigenvals[1])
    end
    map!(p.attributes, [:projs], :nprojs) do projs
        p.nkpts[] == length(projs) ||
            error("Length of projs does not match length of kpoints")
        p.nbands[] == size(projs[1], 1) ||
            error("Length of projs does not match number of bands")
        nprojs = size(projs[1], 2)
        nprojs == length(p.labels[]) ||
            error("Length of projs does not match length of labels")
        return nprojs
    end
    # We use total projectability as the markersize for the scatter plot.
    # Indexed by `band_markersizes[ib][ik]`.
    map!(p.attributes, [:projs, :markersize], :band_markersizes) do projs, markersize
        pj = Wannier.sum_projectability(projs)
        return map(1:p.nbands[]) do ib
            markersize .* [pj[ik][ib] for ik in 1:p.nkpts[]]
        end
    end
    # Get colors for the orbitals from the color scheme
    map!(p.attributes, [:colors, :colorscheme], :orbital_colors) do colors, colorscheme
        colors = isnothing(colors) ? to_colormap(colorscheme) : colors
        # Cycle through the colors if there are more orbitals than colors in the scheme
        return [colors[mod1(ip, length(colors))] for ip in 1:p.nprojs[]]
    end
    # Blend the color of all the orbitals, indexed by `band_colors[ib][ik]`
    map!(p.attributes, [:projs], :band_colors) do projs
        return map(1:p.nbands[]) do ib
            map(1:p.nkpts[]) do ik
                αs = projs[ik][ib, :]
                cs = p.orbital_colors[]
                # Weighted average, such that the order of plotting the
                # orbitals does not matter, and we only need one single
                # scatter plot for each band (instead of nprojs scatters).
                return sum(cs .* αs) / sum(αs)
            end
        end
    end

    # Pass the shared attributes
    kwargs = Dict(k => p[k][] for (k, v) in shared_bandplot_attributes().d)
    # Use a alpha for the color of the bands as background below scatter points
    pb = bandplot!(p, p.kpath[], p.eigenvals[]; kwargs..., linecolor=(:black, 0.3))

    for (ib, bs) in enumerate(pb.bands[])
        scatter!(
            p,
            pb.x[],
            bs;
            markersize=p.band_markersizes[][ib],
            color=p.band_colors[][ib],
        )
    end

    return p
end

function Wannier.get_projbandplot(
    kpath::Wannier.RecipPath, eigenvals::E, projs::P, labels::L;
    show_legend::Bool=true, kwargs...,
) where {
    E<:AbstractVector{<:AbstractVector{<:Real}},
    P<:AbstractVector{<:AbstractMatrix{<:Real}},
    L<:AbstractVector{<:AbstractString},
}
    fig, ax = fig_ax_bandplot(kpath; kwargs...)
    p = projbandplot!(ax, kpath, eigenvals, projs, labels; kwargs...)

    if show_legend
        marker = :circle
        markersize = p.markersize[]
        markers = map(p.orbital_colors[]) do color
            MarkerElement(; color, marker, markersize)
        end
        axislegend(ax, markers, labels)
    end
    return Makie.FigureAxisPlot(fig, ax, p)
end
