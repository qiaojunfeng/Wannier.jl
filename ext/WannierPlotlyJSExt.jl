module WannierPlotlyJSExt
using DocStringExtensions
using Brillouin: KPathInterpolant
using Wannier
using PlotlyJS

"""
    $(SIGNATURES)

Merge consecutive high-symmetry points.

If two high-symmetry kpoints are neighbors, merge them into one,
with label `X|Y`, where `X` and `Y` are the original labels of
the two kpoints, respectively.

# Arguments
- `symm_point_indices`: indices of high-symmetry kpoints, start from 1
- `symm_point_labels`: labels of high-symmetry kpoints
"""
function merge_consecutive_labels(
    symm_point_indices::AbstractVector{<:Integer},
    symm_point_labels::AbstractVector{<:String},
)
    idxs = copy(symm_point_indices)
    labels = copy(symm_point_labels)

    counter = 2
    for i in eachindex(symm_point_indices)[2:end]
        if symm_point_indices[i] == symm_point_indices[i - 1] + 1
            labels[counter - 1] *= "|$(symm_point_labels[i])"
        else
            idxs[counter] = symm_point_indices[i]
            labels[counter] = symm_point_labels[i]
            counter += 1
        end
    end

    return idxs[1:(counter - 1)], labels[1:(counter - 1)]
end

"""
    $(SIGNATURES)

Return a PlotlyJS `Plot` struct for the band structure.

# Arguments
- `x`: Vector for x axis
- `eigenvalues`: band energies, length-`n_kpoints` vector, each element is a
    length-`n_bands` vector

# Keyword Arguments
- `fermi_energy`: Fermi energy, will draw a horizontal line
- `shift_fermi`: shift the Fermi energy to 0
- `symm_point_indices`: indices of high-symmetry kpoints
- `symm_point_labels`: labels of high-symmetry kpoints
- `color`: color of the band, can be
    - a string for the color of all the eigenvalues, e.g., `"red"`
    - a Vector having the same size as `eigenvalues`, to color each individual
        eigenvalue
- `kwargs`: additional keyword arguments for PlotlyJS, the `line` attribute of `scatter`.
    For a full customization, directly manipulate the returned `Plot` struct.
    See [PlotlyJS documentation](http://juliaplots.org/PlotlyJS.jl/stable/syncplots/).
"""
function get_band_plot(
    x::AbstractVector{<:Real},
    eigenvalues::AbstractVector;
    fermi_energy::Union{Nothing,Real}=nothing,
    shift_fermi::Bool=false,
    symm_point_indices::Union{Nothing,AbstractVector}=nothing,
    symm_point_labels::Union{Nothing,AbstractVector}=nothing,
    color="black",
    kwargs...,
)
    nkpts = length(eigenvalues)
    @assert nkpts > 0 "empty eigenvalues"
    nbands = length(eigenvalues[1])

    if symm_point_indices !== nothing
        length(symm_point_indices) == length(symm_point_labels) ||
            error("symm_idx and symm_label must have the same length")
    end
    if shift_fermi && fermi_energy === nothing
        error("shift_fermi is true, but fermi_energy is not given")
    end

    ylabel = "E (eV)"
    # convert to dense matrix for fast indexing, size = nbands * nkpts
    # I use captial letters for matrices
    E = reduce(hcat, eigenvalues)
    if shift_fermi
        E .-= fermi_energy
        ylabel = "E - E_F (eV)"
    end

    if color isa AbstractVector
        uniform_color = false
        # size = nbands * nkpts
        C = reduce(hcat, color)
        cmin, cmax = extrema(C)
        color_plot = Vector(eachrow(C))
        mode = "markers"
        color_kwargs = (; cmin, cmax, colorbar=(;), colorscale="RdBu")
    else
        uniform_color = true
        color_plot = [color for _ in 1:nbands]
        mode = "lines"
    end

    traces = PlotlyJS.AbstractTrace[]
    for (e, c) in zip(eachrow(E), color_plot)
        if uniform_color
            t = PlotlyJS.scatter(; x, y=e, mode, line=(; color=c, kwargs...))
        else
            t = PlotlyJS.scatter(;
                x, y=e, mode, marker=(; color=c, color_kwargs..., kwargs...)
            )
        end
        push!(traces, t)
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
        paper_bgcolor="rgba(0,0,0,0)",  # transparent background
    )

    # for storing infinite lines
    shapes = []

    if symm_point_indices !== nothing
        idxs, labels = merge_consecutive_labels(symm_point_indices, symm_point_labels)
        # add vertial lines for high-symm points to the background
        for i in idxs
            push!(shapes, vline(x[i]; mode="lines", line=attr(; color="black", width=0.2)))
        end
        # labels on x axis
        relayout!(
            layout;
            xaxis=attr(; tickmode="array", tickvals=[x[i] for i in idxs], ticktext=labels),
        )
    end

    if fermi_energy !== nothing
        if shift_fermi
            εF_plot = 0
        else
            εF_plot = fermi_energy
        end
        # add horizontal line for Fermi to the background
        push!(
            shapes,
            hline(εF_plot; mode="lines", line=attr(; dash="dash", color="blue", width=0.2)),
        )
    end

    relayout!(layout; shapes)
    return Plot(traces, layout)
end

function get_band_plot(kpi::KPathInterpolant, eigenvalues::AbstractVector; kwargs...)
    x = Wannier.get_linear_path(kpi)
    symm_point_indices, symm_point_labels = Wannier.get_symm_point_indices_labels(kpi)
    return get_band_plot(x, eigenvalues; symm_point_indices, symm_point_labels, kwargs...)
end

"""
    $(SIGNATURES)

Plot band structure.

# Arguments
- `k`: can be either
    - a Vector of Float64 for x axis value
    - a `KPathInterpolant`
- `eigenvalues`: band energies, length-`n_kpoints` vector, each element is a
    length-`n_bands` vector

# Keyword Arguments
See also the keyword arguments of [`get_band_plot`](@ref).
"""
function Wannier.plot_band(k, eigenvalues::AbstractVector; kwargs...)
    P = get_band_plot(k, eigenvalues; kwargs...)
    return PlotlyJS.plot(P)
end

"""
    $(SIGNATURES)

Compare two band structures.

`eigenvalues_1` in grey, while `eigenvalues_2` in dashed red.

# Arguments
- `k`: can be either
    - a Vector of Float64 for x axis value
    - a `KPathInterpolant`
- `eigenvalues_1`: band energies, see [`plot_band`](@ref)
- `eigenvalues_2`: band energies, see [`plot_band`](@ref)

# Keyword Arguments
See [`plot_band`](@ref)
"""
function Wannier.plot_band_diff(
    k, eigenvalues_1::AbstractArray, eigenvalues_2::AbstractArray; kwargs...
)
    P1 = get_band_plot(k, eigenvalues_1; color="grey", kwargs...)
    # red and slightly thinner
    P2 = get_band_plot(k, eigenvalues_2; color="red", dash="dash", width=0.9, kwargs...)
    addtraces!(P1, P2.data...)
    return plot(P1)
end

end
