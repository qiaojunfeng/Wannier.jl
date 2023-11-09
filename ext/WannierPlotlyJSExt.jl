module WannierPlotlyJSExt
using DocStringExtensions
using Brillouin: KPathInterpolant
using Wannier
using PlotlyJS
using LinearAlgebra
using Printf: @sprintf
using PeriodicTable: elements

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
function Wannier.get_band_plot(
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

function Wannier.get_band_plot(
    kpi::KPathInterpolant, eigenvalues::AbstractVector; kwargs...
)
    x = Wannier.get_linear_path(kpi)
    symm_point_indices, symm_point_labels = Wannier.get_symm_point_indices_labels(kpi)
    return Wannier.get_band_plot(
        x, eigenvalues; symm_point_indices, symm_point_labels, kwargs...
    )
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
See also the keyword arguments of [`Wannier.get_band_plot`](@ref).
"""
function Wannier.plot_band(k, eigenvalues::AbstractVector; kwargs...)
    P = Wannier.get_band_plot(k, eigenvalues; kwargs...)
    return PlotlyJS.plot(P)
end

function Wannier.get_band_diff_plot(
    k, eigenvalues_1::AbstractArray, eigenvalues_2::AbstractArray; kwargs...
)
    P1 = Wannier.get_band_plot(k, eigenvalues_1; color="grey", kwargs...)
    # red and slightly thinner
    P2 = Wannier.get_band_plot(
        k, eigenvalues_2; color="red", dash="dash", width=0.9, kwargs...
    )
    addtraces!(P1, P2.data...)
    return P1
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
    P = Wannier.get_band_diff_plot(k, eigenvalues_1, eigenvalues_2; kwargs...)
    return plot(P)
end

"""A plotly sphere."""
function _sphere(r₀, r; kwargs...)
    N = 32
    u = range(0, 2π, N)
    v = range(0, π, N)
    x = r₀[1] .+ r * cos.(u) * sin.(v)'
    y = r₀[2] .+ r * sin.(u) * sin.(v)'
    z = r₀[3] .+ r * repeat(cos.(v)'; outer=[N, 1])

    return PlotlyJS.surface(; x=x, y=y, z=z, kwargs...)
end

"""
Return a vector of PlotlyJS traces for atoms.

- `lattice`: columnwise lattice vectors
- `atom_positions`: atomic positions in fractional coordinates
- `atom_numbers`: can be either `:H`, `1`, or `"H"`
"""
function _atoms(
    lattice::AbstractMatrix, atom_positions::AbstractVector, atom_numbers::AbstractVector
)
    xyz = []
    radius = []
    color = []
    label = []
    for (pos, atom) in zip(atom_positions, atom_numbers)
        # TODO these are from json file
        # ele = elements[atom]
        # s = ele["symbol"]
        # c = ele["cpkHexColor"]
        # r = ele["radius"] / elements[1]["radius"] / 10  # normalized by Hydrogen radius
        # these from PeriodicTable
        ele = elements[Symbol(atom)]
        s = ele.symbol
        push!(label, s)
        c = ele.cpk_hex
        push!(color, c)
        # https://github.com/JuliaPhysics/PeriodicTable.jl/issues/34
        r = 0.5  # currently no radius in PeriodicTable
        push!(radius, r)
        pos = lattice * pos
        push!(xyz, pos)
    end
    # return xyz, radius, color, label

    # atoms
    # d = Dict(
    #     "mode" => "markers",
    #     "x" => pos[1, :],
    #     "y" => pos[2, :],
    #     "z" => pos[3, :],
    #     "marker" => Dict(
    #         "size" => 5,
    #         "color" => "orange",
    #         # colorscale => "Greens",
    #         # cmin => -20,
    #         # cmax => 50
    #     )
    # )
    # atoms = scatter3d(d)
    atoms = map(atom_positions, xyz, radius, color, label) do frac, cart, r, c, s
        text = @sprintf "%s<br>frac: (%.4f, %.4f, %.4f)</br>cart: (%.4f, %.4f, %.4f)" s frac... cart...
        _sphere(cart, r; colorscale=[[0, c], [1, c]], showscale=false, hovertemplate=text)
    end
    return atoms
end

"""Loop through all the lattice parallelepipede."""
function _lattice_lines(lattice::AbstractMatrix)
    xyz =
        [
            0 0 0
            1 0 0
            1 1 0
            0 1 0
            0 0 0
            0 0 1
            1 0 1
            1 0 0
            1 1 0
            1 1 1
            1 0 1
            1 1 1
            0 1 1
            0 1 0
            0 1 1
            0 0 1
        ]'
    return lattice * xyz
end

"""
Return a vector of PlotlyJS traces for lattice vectors.

- `lattice`: each column is a lattice vector
- `origin`: the overall shift of the structure, in cartesian
"""
function _lattice(lattice::AbstractMatrix, origin::AbstractVector=[0, 0, 0])
    # lattice
    xyz = _lattice_lines(lattice)
    x = xyz[1, :] .+ origin[1]
    y = xyz[2, :] .+ origin[2]
    z = xyz[3, :] .+ origin[3]

    d = Dict(
        "mode" => "lines",
        "x" => x,
        "y" => y,
        "z" => z,
        "line" => Dict("width" => 6, "color" => "black"),
        "hoverinfo" => "none",
    )
    lat = PlotlyJS.scatter3d(d)

    hovers = map(enumerate(eachcol(lattice))) do (i, l)
        @sprintf "a%d<br>(%.4f, %.4f, %.4f)</br>norm: %.4f" i l... norm(l)
    end
    arrow = Dict(
        "x" => lattice[1, :] .+ origin[1],
        "y" => lattice[2, :] .+ origin[2],
        "z" => lattice[3, :] .+ origin[3],
        "u" => 1 .* lattice[1, :],
        "v" => 1 .* lattice[2, :],
        "w" => 1 .* lattice[3, :],
        "anchor" => "tip", # make cone tip be at endpoint
        "sizemode" => "absolute",
        "sizeref" => 0.5,
        # "hoverinfo" => ["text"],
        # "text" => ["a1", "a2", "a3"],
        "hovertemplate" => hovers,
        "colorscale" => [[0, "#FDE74C"], [1, "#FDE74C"]], # color all cones yellow
        "showscale" => false,
    )
    axs = PlotlyJS.cone(arrow)

    return [lat, axs]
end

"""
Return a `PlotlyJS.Plot` for lattice and atoms.

# Arguments
- `lattice`: each column is a lattice vector
- `atom_positions`: each column is an atomic position, in fractional coordinates
- `atom_numbers`: atomic number or symbol of each atom

# Keyword Arguments
- `origin`: the overall shift of the structure, in cartesian
"""
function Wannier.get_crystal_structure_plot(
    lattice::AbstractMatrix,
    atom_positions::AbstractVector,
    atom_numbers::AbstractVector;
    origin::AbstractVector=[0, 0, 0],
)
    traces = _lattice(lattice, origin)

    atoms = _atoms(lattice, atom_positions, atom_numbers)
    append!(traces, atoms)

    ax_attr = attr(;
        tickvals=[],
        zeroline=false,
        showgrid=false,
        showbackground=false,
        title=attr(; text=""),
    )
    layout = Layout(;
        showlegend=false,
        scene=attr(;
            xaxis=ax_attr,
            yaxis=ax_attr,
            zaxis=ax_attr,
            aspectmode="data",
            # camera=attr(;
            #     up=attr(; x=0, z=1, y=0),
            #     center=attr(; x=0, y=0, z=0),
            #     projection=attr(; type="orthographic"),
            # ),
        ),
        # margin=attr(l=0, r=0, b=0, t=0),
        # autosize=false,
        # width=1200, height=1200,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )

    # now get the correct axis range, so that the aspect ratio is correct
    # I need to plot two invisible points representing the bounding box,
    # so that the aspect ratio computed by `aspectmode="data"` is correct.
    # note I use the [:,1] to drop the dim of the returned SMatrix
    center_lattice = 0.5 * sum(lattice; dims=2)[:, 1]
    lattice_extent = _lattice_lines(lattice) .- center_lattice
    atoms_extent = reduce(hcat, [lattice * a for a in atom_positions]) .- center_lattice
    ax_min = minimum(lattice_extent; dims=2)[:, 1]
    ax_max = maximum(lattice_extent; dims=2)[:, 1]
    atoms_min = minimum(atoms_extent; dims=2)[:, 1]
    atoms_max = maximum(atoms_extent; dims=2)[:, 1]
    xyz_min = min.(ax_min, atoms_min)
    xyz_max = max.(ax_max, atoms_max)
    cube_half_width = maximum(abs.([xyz_min..., xyz_max...]))
    # println(center_lattice)
    # println(cube_half_width)
    # I can either put a bounding box
    # o = center_lattice - [cube_half_width, cube_half_width, cube_half_width]
    # a = 2 * [cube_half_width, cube_half_width, cube_half_width]
    # bouding_box = _lattice(Diagonal(a), o)
    # append!(traces, bouding_box)
    # Or just put two points
    x, y, z = zip(center_lattice .- cube_half_width, center_lattice .+ cube_half_width)
    push!(
        traces,
        scatter3d(;
            x=x, y=y, z=z, mode="markers", marker=attr(; size=1, color="rgba(0,0,0,0)")
        ),
    )

    return Plot(traces, layout)
end

"""
Plot crystal structure.
"""
function Wannier.plot_crystal_structure(
    lattice::AbstractMatrix,
    atom_positions::AbstractVector,
    atom_numbers::AbstractVector;
    kwargs...,
)
    P = Wannier.get_crystal_structure_plot(lattice, atom_positions, atom_numbers; kwargs...)
    return PlotlyJS.plot(P)
end

function Wannier.get_dos_plot(
    x::AbstractVector,
    y::AbstractVector;
    fermi_energy=nothing,
    shift_fermi::Bool=false,
    xlabel="Energy (eV)",
    ylabel="DOS (states/eV)",
)
    traces = PlotlyJS.AbstractTrace[]

    t = PlotlyJS.scatter(;
        x,
        y,
        mode="lines",
        # line=(; color=c, kwargs...)
    )
    push!(traces, t)

    layout = PlotlyJS.Layout(;
        showlegend=false,
        xaxis=PlotlyJS.attr(;
            range=[x[1], x[end]],
            title=xlabel,
            zeroline=false,
            showgrid=true,
            showbackground=false,
            # gridcolor = KLINE_COL[],
            ticks="outside",
            showline=true,
            mirror=true,
            linecolor="black", # show axis boundary
        ),
        yaxis=PlotlyJS.attr(;
            range=[minimum(y) - 0.5, maximum(y) + 0.5],
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

    shapes = []
    # add horizontal line for y = 0
    push!(
        shapes,
        PlotlyJS.hline(
            0;
            mode="lines",
            line=PlotlyJS.attr(; dash="dash", color="grey", width=0.2),
        ),
    )
    if fermi_energy !== nothing
        if shift_fermi
            x .-= fermi_energy
            ylabel = "E - E_F (eV)"
            εF_plot = 0
        else
            εF_plot = fermi_energy
        end
        # add vertical line for Fermi to the background
        push!(
            shapes,
            PlotlyJS.vline(
                εF_plot;
                mode="lines",
                line=PlotlyJS.attr(; dash="dash", color="blue", width=0.2),
            ),
        )
    end
    PlotlyJS.relayout!(layout; shapes)

    return PlotlyJS.Plot(traces, layout)
end

"""
Plot density of states.
"""
function Wannier.plot_dos(x::AbstractVector, y::AbstractVector; kwargs...)
    P = Wannier.get_dos_plot(x, y; kwargs...)
    return PlotlyJS.plot(P)
end

end
