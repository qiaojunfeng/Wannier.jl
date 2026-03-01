module WannierMakieExt

using Makie

using Wannier
import Wannier: bandplot, bandplot!

@recipe BandPlot begin
    "Fermi energy"
    fermi_energy = nothing

    "Whether to shift the bands such that the Fermi energy is at zero"
    shift_fermi = false

    "Legend label for the bands"
    label = nothing

    "Whether to show the legend"
    show_legend = false

    "Line color"
    linecolor = @inherit linecolor :black

    "Line style"
    linestyle = @inherit linestyle :solid

    "Line width"
    linewidth = @inherit linewidth 1

    "Line color for vertical lines at high-symmetry kpoints"
    linecolor_vline = :grey

    "Line style for vertical lines"
    linestyle_vline = :solid

    "Vertical line width"
    linewidth_vline = 0.8

    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(
    p::BandPlot{Tuple{X,T,L,E}}
) where {
    X<:AbstractVector{<:Real},T<:AbstractVector{<:Integer},
    L<:AbstractVector{<:AbstractString},E<:AbstractVector{<:AbstractVector{<:Real}}
}
    # X: cumulative distance along kpath
    x = p[1][]
    # T: indices of high-symmetry kpoints
    xt_idxs = p[2][]
    # L: labels of high-symmetry kpoints
    xt_labels = p[3][]
    # E: eigenvalues of the bands
    eigenvals = p[4][]

    nkpts = length(eigenvals)
    nkpts == length(x) || error("Length of eigenvals does not match length of kpoints")
    nbands = length(eigenvals[1])

    map!(
        p.attributes, [:fermi_energy, :shift_fermi], :bands
    ) do fermi_energy, shift_fermi
        # The input eigenvals is a vector of vectors, where the outer vector is
        # over k-points and the inner vector is over bands. We need to reshape
        # it to a vector of length nbands with inner vector of length nkpoints
        # for plotting.
        bands = map(1:nbands) do ib
            [eigenvals[ik][ib] for ik in 1:nkpts]
        end
        if shift_fermi && !isnothing(fermi_energy)
            return bands .- fermi_energy
        else
            return bands
        end
    end

    map!(p.attributes, [:shift_fermi], :ylabel) do shift_fermi
        if shift_fermi
            return rich(rich("E - E"; font=:italic), subscript("F"), " (eV)")
        else
            return "Energy (eV)"
        end
    end

    xts = x[xt_idxs]
    # xticks!(parent_scene; xtickrange=xts, xticklabels=xt_labels)
    vlines!(p, xts; color=p.linecolor_vline[], linewidth=p.linewidth_vline[], linestyle=p.linestyle_vline[])

    for (i, bs) in enumerate(p.bands[])
        kwargs = (; color=p.linecolor[], linestyle=p.linestyle[], linewidth=p.linewidth[])
        # Only show label for the first band, to avoid cluttering the legend
        (i == 1) && (kwargs = (; kwargs..., label=p.label[]))
        lines!(p, x, bs; kwargs...)
    end

    # This is a bit hacky: the recipe only works with the Plot, but we want to
    # set the ylabel of the Axis.

    # 1. Find the Axis containing this plot
    # p.parent is usually the Scene; we check the Scene's parent
    parent_scene = p.parent
    parent_ax = parent_scene.parent
    @info "ab" parent_scene parent_ax

    # 2. Safety check: Ensure we are actually in an Axis (not just a raw Scene)
    if parent_ax isa Axis
        # 3. Only set if the current label is empty to avoid overwriting the user
        if isempty(parent_ax.ylabel[])
            parent_ax.ylabel = p.ylabel[]
        end
    end

    p.show_legend[] && legend!(p)

    return p
end

function Makie.convert_arguments(
    ::Type{<:BandPlot}, kpath::Wannier.RecipPath, eigenvals::E
) where {E<:AbstractVector{<:AbstractVector{<:Real}}}
    x = Wannier.get_linear_path(kpath)
    labs = Wannier.symm_label_to_unicode(kpath.labels)
    idxs, labs = Wannier.merge_symm_labels(kpath.indices, labs)
    return (x, idxs, labs, eigenvals)
end

end # module
