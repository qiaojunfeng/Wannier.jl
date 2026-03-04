import Wannier: projbandplot, projbandplot!, get_projbandplot

@recipe ProjBandPlot (kpath, eigenvals, projs, labels) begin
    shared_bandplot_attributes()...

    "Marker size for the scatter points representing the projectabilities"
    markersize = @inherit markersize 10

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
    # Projectability for each eigenstate onto each orbital, accessed by P[ip][ib][ik]
    map!(p.attributes, [:projs], :band_projs) do ps
        return map(1:p.nprojs[]) do ip
            map(1:p.nbands[]) do ib
                map(1:p.nkpts[]) do ik
                    ps[ik][ib, ip]
                end
            end
        end
    end
    # We use total projectability as the markersize for the scatter plot.
    # Accessed by P[ib][ik].
    map!(p.attributes, [:band_projs], :band_proj_tot) do proj
        return sum(proj)
    end
    # Get colors for the orbitals from the color scheme
    map!(p.attributes, [:colorscheme], :orbital_colors) do colorscheme
        colors = to_colormap(colorscheme)
        # Cycle through the colors if there are more orbitals than colors in the scheme
        return [colors[mod1(ip, length(colors))] for ip in 1:p.nprojs[]]
    end

    # Pass the shared attributes
    kwargs = Dict(k => p[k][] for (k, v) in shared_bandplot_attributes().d)
    # Use a alpha for the color of the bands as background below scatter points
    pb = bandplot!(p, p.kpath[], p.eigenvals[]; kwargs..., linecolor=(:black, 0.3))

    for (ip, ps) in enumerate(p.band_projs[])
        # Size scale of markers
        s = p.markersize[]
        # Color for this orbital
        c = RGBf(p.orbital_colors[][ip])
        for (ib, bs) in enumerate(pb.bands[])
            # Only show label for the first band, to avoid cluttering the legend
            kwargs = ib == 1 ?  (; label=p.labels[][ip]) : (;)
            scatter!(
                p,
                pb.x[],
                bs;
                markersize=s * p.band_proj_tot[][ib],
                # Color with some transparency
                color=[RGBAf(c, α) for α in ps[ib]],
                kwargs...,
            )
        end
    end

    return p
end

function Wannier.get_projbandplot(
    kpath::Wannier.RecipPath, eigenvals::E, projs::P, labels::L; kwargs...
) where {
    E<:AbstractVector{<:AbstractVector{<:Real}},
    P<:AbstractVector{<:AbstractMatrix{<:Real}},
    L<:AbstractVector{<:AbstractString},
}
    fig, ax = fig_ax_bandplot(kpath; kwargs...)
    p = projbandplot!(ax, kpath, eigenvals, projs, labels; kwargs...)

    if length(labels) > 1
        scts = filter(p.plots) do p
            l = get(p.attributes, :label, nothing)
            (p isa Scatter) && !isnothing(l)
        end
        scts = map(scts) do s
            # Remove transparency for the legend
            color = RGBAf(RGBf(s.color[][1]), 1.0)
            s => (; color)
        end
        axislegend(ax, scts, labels)
    end
    return Makie.FigureAxisPlot(fig, ax, p)
end
