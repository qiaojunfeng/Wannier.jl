# Must prepend a dot due to Requires.jl
using .PlotlyJS

export plot_band

"""
    _merge_consecutive_labels(
        symm_idx::AbstractArray{Int}, symm_label::AbstractVector{String}
    )

Merge consecutive high-symmetry points.

If two high-symmetry kpoints are neighbors, merge them into one,
with label `X|Y`, where `X` and `Y` are the original labels of
the two kpoints, respectively.

# Arguments
- `symm_idx`: indices of high-symmetry kpoints, start from 1
- `symm_label`: labels of high-symmetry kpoints
"""
function _merge_consecutive_labels(
    symm_idx::AbstractArray{Int}, symm_label::AbstractVector{String}
)
    idx = copy(symm_idx)
    labels = copy(symm_label)

    counter = 2
    for i in axes(symm_idx, 1)[2:end]
        if symm_idx[i] == symm_idx[i - 1] + 1
            labels[counter - 1] *= "|$(symm_label[i])"
        else
            idx[counter] = symm_idx[i]
            labels[counter] = symm_label[i]
            counter += 1
        end
    end

    return idx[1:(counter - 1)], labels[1:(counter - 1)]
end

"""
    _get_band_plot(
        x::AbstractVector{T},
        E::AbstractArray{T};
        fermi_energy::Union{Nothing,T}=nothing,
        shift_fermi::Bool=false,
        symm_idx::Union{Nothing,AbstractArray{Int}}=nothing,
        symm_label::Union{Nothing,AbstractVector{String}}=nothing,
        color="black",
        kwargs...,
    ) where {T<:Real}

Return a PlotlyJS `Plot` struct for the band structure.

# Arguments
- `x`: 1D array for x axis
- `E`: band energies, 1D of length `n_kpts`, or 2D array of size `n_bands * n_kpts`

# Keyword Arguments
- `fermi_energy`: Fermi energy, will draw a horizontal line
- `shift_fermi`: shift the Fermi energy to 0
- `symm_idx`: indices of high-symmetry kpoints
- `symm_label`: labels of high-symmetry kpoints
- `color`: color of the band
- `kwargs`: additional keyword arguments for PlotlyJS, the `line` attribute of `scatter`.
    For a full customization, directly manipulate the returned `Plot` struct.
    See [PlotlyJS documentation](http://juliaplots.org/PlotlyJS.jl/stable/syncplots/).
"""
function _get_band_plot(
    x::AbstractVector{T},
    E::AbstractArray{T};
    fermi_energy::Union{Nothing,T}=nothing,
    shift_fermi::Bool=false,
    symm_idx::Union{Nothing,AbstractArray{Int}}=nothing,
    symm_label::Union{Nothing,AbstractVector{String}}=nothing,
    color="black",
    kwargs...,
) where {T<:Number}
    ndims(E) <= 2 || error("E must be a 1D or 2D array")
    if symm_idx !== nothing
        length(symm_idx) == length(symm_label) ||
            error("symm_idx and symm_label must have the same length")
    end
    if shift_fermi && fermi_energy === nothing
        error("shift_fermi is true, but fermi_energy is not given")
    end

    ylabel = "E (eV)"
    E_plot = E
    if shift_fermi
        E_plot = E - fermi_energy
        ylabel = "E - E_F (eV)"
    end

    if ndims(E) == 2
        traces = PlotlyJS.AbstractTrace[]
        for e in eachrow(E_plot)
            push!(traces, scatter(; x=x, y=e, line=attr(; color=color, kwargs...)))
        end
    else
        traces = scatter(; x=x, y=E_plot, line=attr(; color=color, kwargs...))
    end

    layout = Layout(;
        showlegend=false,
        xaxis=attr(;
            range=[x[1], x[end]],
            zeroline=false,
            showgrid=true,
            showbackground=false,
            # gridcolor = KLINE_COL[],
            ticks="outside",
            showline=true,
            mirror=true,
            linecolor="black", # show axis boundary
        ),
        yaxis=attr(;
            range=[minimum(E) - 0.5, maximum(E) + 0.5],
            title=ylabel,
            zeroline=false,
            showgrid=false,
            showbackground=false,
            ticks="outside",
            showline=true,
            mirror="all",
            linecolor="black", # show axis boundaries on all subplots
        ),
        hovermode="closest",
        autosize=true,
        # width = 480, height = 480,
        # #margin=attr(l=50, r=5, b=15, t=10),
        plot_bgcolor="rgba(0,0,0,0)",  # transparent background
        paper_bgcolor="rgba(0,0,0,0)",  # transparent backgroun
    )

    # for storing infinite lines
    shapes = []

    if symm_idx !== nothing
        idx, label = _merge_consecutive_labels(symm_idx, symm_label)
        # add vertial lines for high-symm points to the background
        for i in idx
            push!(shapes, vline(x[i]; mode="lines", line=attr(; color="black", width=0.2)))
        end
        # labels on x axis
        relayout!(
            layout;
            xaxis=attr(; tickmode="array", tickvals=[x[i] for i in idx], ticktext=label),
        )
    end

    if fermi_energy !== nothing
        # add horizontal line for Fermi to the background
        push!(
            shapes,
            hline(
                fermi_energy;
                mode="lines",
                line=attr(; dash="dash", color="blue", width=0.2),
            ),
        )
    end

    relayout!(layout; shapes=shapes)

    return Plot(traces, layout)
end

"""
    plot_band(x::AbstractVector, E::AbstractArray; kwargs...)

Plot band structure.

# Arguments
- `x`: 1D array for x axis
- `E`: band energies, 1D of length `n_kpts`, or 2D array of size `n_bands * n_kpts`

# Keyword Arguments
See also the keyword arguments of [`_get_band_plot`](@ref).
"""
function plot_band(x::AbstractVector, E::AbstractArray; kwargs...)
    P = _get_band_plot(x, E; kwargs...)
    return plot(P)
end

"""
    plot_band(kpi::KPathInterpolant, E::AbstractArray; kwargs...)

Plot band structure.

# Arguments
- `kpi`: KPathInterpolant
- `E`: band energies, 1D of length `n_kpts`, or 2D array of size `n_bands * n_kpts`

# Keyword Arguments
See also the keyword arguments of [`_get_band_plot`](@ref).
"""
function plot_band(kpi::KPathInterpolant, E::AbstractArray; kwargs...)
    x = get_x(kpi)
    symm_idx, symm_label = _get_symm_idx_label(kpi)
    return plot_band(x, E; symm_idx=symm_idx, symm_label=symm_label, kwargs...)
end
