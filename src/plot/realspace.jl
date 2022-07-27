using LinearAlgebra
using GeometryBasics
using Meshing
using PlotlyJS
using GLMakie
# using PeriodicTable
using JSON
# using Colors
using StatsBase
using Wannier

# I haven't include this file into the package, because this introduce lots of dependencies.
# TODO Once I sorted out how to incoporate plotting functions into the package without significantly
# slowing down the precompilation, I will include this file into the package.

# unfortunately, the PeriodicTable.jl does not provide atomic radius info,
# we can use this json from this repo
# https://github.com/AlexGustafsson/molecular-data/blob/master/json/elements.json
# TODO once there is a julia package providing atomic radius and color, switch to that one
# You need to first run
#   wget https://raw.githubusercontent.com/AlexGustafsson/molecular-data/master/json/elements.json
elements = JSON.parsefile(joinpath(@__DIR__, "elements.json"))

function _plotly_mesh3d(
    origin::AbstractVector{T},
    lattice::AbstractMatrix{T},
    W::AbstractArray{T,3},
    iso::T,
    color,
) where {T<:Real}
    algo = MarchingCubes(; iso=iso, insidepositive=true)
    # use marching tetrahedra with iso at 100
    # algo = MarchingTetrahedra(iso=100, insidepositive=true)
    # use Naive Surface Nets with iso at 100
    # algo = NaiveSurfaceNets(iso=100, insidepositive=true)

    # generate the mesh using marching cubes
    # mc = Mesh(cube.W, algo)
    # save the file as a PLY file (change extension to save as STL, OBJ, OFF)
    # save("mc.ply", mc)

    # call isosurface to get a vector of points and vector of faces indexing to the points
    widths = Point{3,Float64}([norm(lattice[:, i]) for i in 1:3])
    # println(origin, widths)
    vertices, faces = Meshing.isosurface(
        W, algo, Point{3,Float32}, TriangleFace{Int}; origin=origin, widths=widths
    )

    n_vert = length(vertices)
    xyz = zeros(Float64, 3, n_vert)

    for i in 1:n_vert
        xyz[:, i] = vertices[i]
    end

    n_face = length(faces)
    ijk = zeros(Int, 3, n_face)
    for i in 1:n_face
        # plotly starts from 0
        ijk[:, i] = faces[i] .- 1
    end

    x = xyz[1, :]
    y = xyz[2, :]
    z = xyz[3, :]
    i = ijk[1, :]
    j = ijk[2, :]
    k = ijk[3, :]

    t = PlotlyJS.mesh3d(; x=x, y=y, z=z, i=i, j=j, k=k, color=color, opacity=0.6)
    return t
end

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
Plotly lattice and atoms

lattice: each column is a lattice vector
origin: the overall shift of the structure
atom_positions: each column is an atomic position
atom_numbers: atomic number of each atom
"""
function _plotly_structure(
    origin::AbstractVector,
    lattice::AbstractMatrix,
    atom_positions::AbstractMatrix,
    atom_numbers::AbstractVector,
)
    # trace to loop through all the lattice parallelepipede
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
    xyz = lattice * xyz
    x = xyz[1, :] .+ origin[1]
    y = xyz[2, :] .+ origin[2]
    z = xyz[3, :] .+ origin[3]

    d = Dict(
        "mode" => "lines",
        "x" => x,
        "y" => y,
        "z" => z,
        "line" => Dict("width" => 6, "color" => "black"),
    )
    lat = PlotlyJS.scatter3d(d)

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
        "hoverinfo" => ["text"],
        "text" => ["a1", "a2", "a3"],
        "colorscale" => [[0, "#FDE74C"], [1, "#FDE74C"]], # color all cones yellow
        "showscale" => false,
    )
    axs = PlotlyJS.cone(arrow)

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
    atoms = []
    for i in 1:size(atom_positions, 2)
        # TODO find a package to get atom info instead of json
        ele = elements[atom_numbers[i]]
        s = ele["symbol"]
        c = ele["cpkHexColor"]
        r = ele["radius"] / elements[1]["radius"] / 10  # normalized by Hydrogen radius
        colorscale = [[0, c], [1, c]]
        push!(
            atoms,
            _sphere(
                atom_positions[:, i], r; colorscale=colorscale, text=s, showscale=false
            ),
        )
    end

    return [lat, axs, atoms...]
end

function _guess_isolevel(data)
    h = fit(Histogram, vec(data); nbins=100)
    percent = cumsum(h.weights) / length(data)
    i = findfirst(percent .>= 0.97)
    # since length(h.edges) = length(h.weights) + 1
    return h.edges[1][i + 1]
end

"""
Plot cube file with Plotly.

data: volumetric data in 3D
"""
function plot_cube_plotly(
    origin::AbstractVector{T},
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{U},
    data::AbstractArray{T,3};
    iso::Union{T,Nothing}=nothing,
) where {T<:Real,U<:Integer}
    structure = _plotly_structure(origin, lattice, atom_positions, atom_numbers)

    if isnothing(iso)
        iso = _guess_isolevel(data)
    end
    a, b = minimum(data), maximum(data)
    if iso >= a
        s1 = _plotly_mesh3d(origin, lattice, data, iso, "#C3423F")
    end
    if iso <= b
        s2 = _plotly_mesh3d(origin, lattice, data, -iso, "#5BC0EB")
    end

    return PlotlyJS.plot([structure..., s1, s2])
end

function plot_cube_makie(
    origin::AbstractVector{T},
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{U},
    data::AbstractArray{T,3};
    iso::Union{T,Nothing}=nothing,
) where {T<:Real,U<:Integer}
    scene = Makie.Scene()
    cam3d!(scene)

    if isnothing(iso)
        iso = _guess_isolevel(data)
    end
    a, b = minimum(data), maximum(data)

    if iso >= a
        algo = MarchingCubes(; iso=iso, insidepositive=true)
        mc = GeometryBasics.Mesh(data, algo)
        Makie.mesh!(scene, mc; color="#C3423F", shading=true)#color=[norm(v) for v in coordinates(mc)])
    end
    if iso <= b
        algo = MarchingCubes(; iso=-iso, insidepositive=true)
        mc = GeometryBasics.Mesh(data, algo)
        Makie.mesh!(scene, mc; color="#5BC0EB", shading=true)#color=[norm(v) for v in coordinates(mc)])
    end

    return scene
end

cube = read_cube(expanduser("~/BaTiO3_00002.cube"))
plot_cube_plotly(cube.origin, cube.lattice, cube.atom_positions, cube.atom_numbers, cube.W)
# plot_cube_makie(cube.origin, cube.lattice, cube.atom_positions, cube.atom_numbers, cube.W)
