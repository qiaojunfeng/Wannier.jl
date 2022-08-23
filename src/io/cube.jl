using Printf: @printf
using Dates: now
using LazyGrids: ndgrid

# Cube format
# Specification from http://paulbourke.net/dataformats/cube/

"""
Read cube file.

All outputs in cartesian coordinates, bohr unit.
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

    atom_positions = zeros(Float64, 3, n_atoms)
    atom_numbers = zeros(Int, n_atoms)
    for i in 1:n_atoms
        line = split(strip(readline(io)))
        atom_numbers[i] = parse(Int, line[1])
        charge = parse(Float64, line[2])
        # cartesian coordinates, Bohr unit
        atom_positions[:, i] = parse.(Float64, line[3:5])
    end

    n_x, n_y, n_z = n_voxels
    X = 0:(n_x - 1)
    Y = 0:(n_y - 1)
    Z = 0:(n_z - 1)
    # fractional w.r.t. span_vectors
    Xg, Yg, Zg = ndgrid(X, Y, Z)

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
    return (
        atom_positions=atom_positions,
        atom_numbers=atom_numbers,
        origin=origin,
        span_vectors=span_vectors,
        X=Xg,
        Y=Yg,
        Z=Zg,
        W=W,
    )
end

"""
Write cube file.

atom_positions: 3 x n_atoms, Å, fractional coordinates
atom_numbers: n_atoms, atomic numbers
origin: 3, Å, origin of the grid
span_vectors: 3 x 3, Å, each column is a spanning vector
W: nx x ny x nz, volumetric data
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
        ax = span_vectors[:, i] ./ Bohr
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
