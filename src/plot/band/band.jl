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

Return a PlotlyJS `Plot` struct for the band structure.

# Arguments
- `x`: 1D array for x axis
- `E`: band energies, 1D of length `n_kpts`, or 2D array of size `n_bands * n_kpts`

# Keyword Arguments
"""
function _get_band_plot(
    x::AbstractVector{T},
    E::AbstractArray{T};
    fermi_energy::Union{Nothing,T}=nothing,
    shift_fermi::Bool=false,
    symm_idx::Union{Nothing,AbstractArray{Int}}=nothing,
    symm_label::Union{Nothing,AbstractVector{String}}=nothing,
    kwargs...,
) where {T<:Real}
    ndims(E) <= 2 || error("E must be a 1D or 2D array")
    if symm_idx !== nothing
        length(symm_idx) == length(symm_label) ||
            error("symm_idx and symm_label must have the same length")
    end
    if shift_fermi && fermi_energy === nothing
        error("shift_fermi is true, but fermi_energy is not given")
    end

    E_plot = E
    if ndims(E) == 2
        E_plot = E'
    end

    ylabel = "E (eV)"
    if shift_fermi
        E_plot -= fermi_energy
        ylabel = "E - E_F (eV)"
    end

    traces = scatter(; x=x, y=E_plot, kwargs...)

    layout = Layout(;
        # showlegend=false,
        xaxis=attr(;
            range=[x[1], x[end]],
            # zeroline=false,
            # showgrid=true, showbackground=false,
            # gridcolor = KLINE_COL[], ticks="outside",
            # showline=true, mirror=true, linecolor="black", # show axis boundary
        ),
        yaxis=attr(;
            range=[minimum(E) - 0.5, maximum(E) + 0.5],
            title=ylabel,
            # zeroline=false,
            # showgrid=false, showbackground=false,
            # ticks="outside",
            # showline=true, mirror="all", linecolor="black", # show axis boundaries on all subplots
        ),
        # hovermode = "closest",
        # autosize = true,
        # width = 480, height = 480,
        # #margin=attr(l=50, r=5, b=15, t=10),
        # plot_bgcolor=TRANSPARENT_COL[], paper_bgcolor=TRANSPARENT_COL[],
    )

    if fermi_energy !== nothing
        # add vertical line for Fermi to the background
        prepend!(
            traces,
            scatter(;
                x=[x[1], x[end]],
                y=[fermi_energy, fermi_energy],
                line=attr(; dash="dash", color="blue", width=0.2),
            ),
        )
    end

    if symm_idx !== nothing
        idx, label = _merge_consecutive_labels(symm_idx, symm_label)
        # add vertial lines for high-symm points to the background
        for i in idx
            prepend!(
                traces,
                scatter(;
                    x=[x[i], x[i]],
                    y=[minimum(E), maximum(E)],
                    line=attr(; color="black", width=0.2, fill="toself"),
                ),
            )
        end
        # labels on x axis
        relayout!(
            layout;
            xaxis=attr(; tickmode="array", tickvals=[x[i] for i in idx], ticktext=label),
        )
    end

    return Plot(traces, layout)
end

function plot_band(x::AbstractVector{T}, E::AbstractArray{T}; kwargs...) where {T<:Real}
    P = _get_band_plot(x, E, kwargs...)
    return plot(P)
end

function plot_band(kpi::KPathInterpolant, E::AbstractArray{T}; kwargs...) where {T<:Real}
    x = get_x(kpi)
    return plot_band(x, E; kwargs...)
end
