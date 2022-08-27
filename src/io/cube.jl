using Printf: @printf
using Dates: now
using LazyGrids: ndgrid

# Cube format
# Specification from http://paulbourke.net/dataformats/cube/

export read_cube, write_cube

"""
    read_cube(filename::AbstractString)

Read `cube` file.

!!! note

    By default, `cube` use Bohr unit, here all returns are in Cartesian coordinates, Å unit.
"""
function read_cube(filename::AbstractString)
    @info "Reading cube file: " filename
    io = open(filename)

    # header
    header = readline(io; keep=true)
    header *= readline(io; keep=true)
    print(header)

    line = split(strip(readline(io)))
    n_atoms = parse(Int, line[1])
    origin = parse.(Float64, line[2:4])
    # to Å unit
    origin .*= Bohr

    # number of voxels per spanning vector
    n_voxels = zeros(Int, 3)
    span_vectors = zeros(Float64, 3, 3)
    for i in 1:3
        line = split(strip(readline(io)))
        n_v = parse(Int, line[1])
        n_voxels[i] = n_v
        # bohr unit
        span_vectors[:, i] = parse.(Float64, line[2:4])
    end
    # to Å unit
    span_vectors .*= Bohr

    atom_positions = zeros(Float64, 3, n_atoms)
    atom_numbers = zeros(Int, n_atoms)
    for i in 1:n_atoms
        line = split(strip(readline(io)))
        atom_numbers[i] = parse(Int, line[1])
        charge = parse(Float64, line[2])
        # cartesian coordinates, Bohr unit
        atom_positions[:, i] = parse.(Float64, line[3:5])
    end
    # to Å unit
    atom_positions .*= Bohr

    n_x, n_y, n_z = n_voxels
    # fractional w.r.t. span_vectors
    X = range(0, 1, n_x)
    Y = range(0, 1, n_y)
    Z = range(0, 1, n_z)
    Xg, Yg, Zg = ndgrid(X, Y, Z)
    # fractional w.r.t. span_vectors
    O = inv(span_vectors) * origin
    X = range(0, 1, n_x) .- O[1]
    Y = range(0, 1, n_y) .- O[2]
    Z = range(0, 1, n_z) .- O[3]
    Xg, Yg, Zg = ndgrid(X, Y, Z)
    rgrid = RGrid(span_vectors, Xg, Yg, Zg)

    W = zeros(Float64, n_x, n_y, n_z)
    # 6 columns per line
    d, r = divrem(n_z, 6)
    if r > 0
        nline = d + 1
    else
        nline = d
    end
    for ix in 1:n_x
        for iy in 1:n_y
            iz = 1
            for _ in 1:nline
                line = split(strip(readline(io)))
                n = length(line)
                # skip empty line
                if n == 0
                    line = split(strip(readline(io)))
                    n = length(line)
                end
                W[ix, iy, iz:(iz + n - 1)] = parse.(Float64, line)
                iz += n
            end
        end
    end

    close(io)
    return (atom_positions=atom_positions, atom_numbers=atom_numbers, rgrid=rgrid, W=W)
end

"""
    write_cube(filename, filename, atom_positions, atom_numbers, origin, span_vectors, W)

Write `cube` file.

# Arguments
- `atom_positions`: `3 * n_atoms`, Å, cartesian coordinates
- `atom_numbers`: `n_atoms`, atomic numbers
- `origin`: `3`, Å, origin of the grid
- `span_vectors`: `3 * 3`, Å, each column is a spanning vector
- `W`: `nx * ny * nz`, volumetric data

See also [`write_cube(filename, lattice, atom_positions, atom_numbers, wf_center, rgrid, W; radius=4.0)`](@ref write_cube(filename, lattice, atom_positions, atom_numbers, wf_center, rgrid, W; radius=4.0))
"""
function write_cube(
    filename::AbstractString,
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{Int},
    origin::AbstractVector{T},
    span_vectors::AbstractMatrix{T},
    W::AbstractArray{T,3},
) where {T<:Real}
    n_atoms = length(atom_numbers)
    size(atom_positions, 2) == n_atoms || error("incompatible n_atoms")
    size(span_vectors) == (3, 3) || error("incompatible span_vectors")
    length(origin) == 3 || error("origin must be 3-vector")

    @info "Writing cube file: " filename
    io = open(filename, "w")

    # header
    @printf(io, "Created by Wannier.jl %s\n", string(now()))
    @printf(io, "outer loop: x, middle loop: y, inner loop: z\n")

    # to Bohr
    origin_bohr = origin ./ Bohr
    @printf(io, "%d %12.6f %12.6f %12.6f\n", n_atoms, origin_bohr...)

    n_xyz = size(W)
    for i in 1:3
        # number of voxels
        n_v = n_xyz[i]
        ax = span_vectors[:, i] ./ (n_v * Bohr)
        @printf(io, "%d %12.6f %12.6f %12.6f\n", n_v, ax...)
    end

    for i in 1:n_atoms
        n = atom_numbers[i]
        charge = 1.0
        pos = atom_positions[:, i] ./ Bohr
        @printf(io, "%d %12.6f %12.6f %12.6f %12.6f\n", n, charge, pos...)
    end

    for ix in 1:n_xyz[1]
        for iy in 1:n_xyz[2]
            for iz in 1:n_xyz[3]
                @printf(io, "%12.6g ", W[ix, iy, iz])

                if (iz % 6 == 0)
                    @printf(io, "\n")
                end
            end
            @printf(io, "\n")
        end
    end

    close(io)
    return nothing
end

"""
    write_cube(filename, lattice, atom_positions, atom_numbers, wf_center, rgrid, W; radius=4.0)

Write `cube` file for WF.

# Arguments
- `lattice`: each column is a lattice vector, Å
- `atom_positions`: `3 * n_atoms`, fractional coordinates
- `atom_numbers`: `n_atoms`, atomic numbers
- `wf_centers`: `3`, fractional coordinates w.r.t. lattice
- `rgrid`: `RGrid`
- `W`: `nx * ny * nz`, volumetric data

# Keyword arguments
- `radius`: Å. Periodic replica of atoms are written only for the region
    within `radius` Å from `wf_center`.

See also [`write_cube(filename, filename, atom_positions, atom_numbers, origin, span_vectors, W)`]
(@ref write_cube(filename, filename, atom_positions, atom_numbers, origin, span_vectors, W))
"""
function write_cube(
    filename::AbstractString,
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{Int},
    wf_centers::AbstractVector{T},
    rgrid::RGrid,
    W::AbstractArray{T,3};
    radius::T=4.0,
) where {T<:Real}
    O = origin(rgrid)
    spanvec = span_vectors(rgrid)
    # move atom closer to WF center
    n_atoms = size(atom_positions, 2)
    atom_c = sum(atom_positions; dims=2) / n_atoms
    # translations for atoms, and to cartesian
    t = round.(wf_centers .- atom_c)
    positions = lattice * (atom_positions .+ t)
    # if atoms from neighbouring cells are close to wfc, include them
    wfc = lattice * wf_centers  # cartesian
    all_pos = [Vector(c) for c in eachcol(positions)]
    all_atoms = [c for c in atom_numbers]
    for i in -3:3
        for j in -3:3
            for k in -3:3
                pos = positions .+ lattice * [i, j, k]
                for n in 1:n_atoms
                    if norm(pos[:, n] .- wfc) <= radius
                        push!(all_atoms, atom_numbers[n])
                        push!(all_pos, pos[:, n])
                    end
                end
            end
        end
    end
    all_pos = reduce(hcat, all_pos)
    # fake atom for WF center
    # all_atoms = push!(all_atoms, 100)
    # all_pos = hcat(all_pos, wfc)

    return write_cube(filename, all_pos, all_atoms, O, spanvec, W)
end
