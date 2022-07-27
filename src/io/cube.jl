using Printf: @printf
using Dates: now

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

    # number of voxels per lattice vector
    n_voxels = zeros(Int, 3)
    lattice = zeros(Float64, 3, 3)
    for i in 1:3
        line = split(strip(readline(io)))
        n_v = parse(Int, line[1])
        n_voxels[i] = n_v
        # bohr unit
        lattice[:, i] = (n_v - 1) * parse.(Float64, line[2:4])
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
    X = 1:n_x
    Y = 1:n_y
    Z = 1:n_z
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
        lattice=lattice,
        X=X,
        Y=Y,
        Z=Z,
        W=W,
    )
end

"""
Write cube file.

atom_positions: 3 x n_atoms, fractional coordinates
atom_numbers: n_atoms, atomic numbers
lattice: 3 x 3, Å, each column is a lattice vector
n_voxels: 3, number of voxels along three lattice vectors
X: fractional coordinates of W along x
Y: fractional coordinates of W along y
Z: fractional coordinates of W along z
W: nx x ny x nz, volumetric data
"""
function write_cube(
    filename::AbstractString,
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{Int},
    lattice::AbstractMatrix{T},
    X::AbstractVector{T},
    Y::AbstractVector{T},
    Z::AbstractVector{T},
    W::AbstractArray{T,3},
) where {T<:Real}
    n_atoms = length(atom_numbers)
    size(atom_positions, 2) != n_atoms && error("incompatible n_atoms")
    size(lattice) != (3, 3) && error("incompatible lattice")

    n_x, n_y, n_z = size(W)
    length(X) != n_x && error("incompatible n_gx")
    length(Y) != n_y && error("incompatible n_gy")
    length(Z) != n_z && error("incompatible n_gz")

    @info "Writing cube file: " filename
    io = open(filename, "w")

    # to Bohr
    lattice_bohr = lattice ./ Bohr

    # header
    @printf(io, "Created by Wannier.jl %s\n", string(now()))
    @printf(io, "outer loop: x, middle loop: y, inner loop: z\n")

    origin = lattice_bohr * [X[1], Y[1], Z[1]]
    @printf(io, "%d %12.6f %12.6f %12.6f\n", n_atoms, origin...)

    # number of voxels per lattice vector
    n_voxels = size(W)
    for i in 1:3
        n_v = n_voxels[i]
        ax = lattice_bohr[:, i] / n_v
        @printf(io, "%d %12.6f %12.6f %12.6f\n", n_v, ax...)
    end

    for i in 1:n_atoms
        n = atom_numbers[i]
        charge = 1.0
        # to Bohr
        pos = lattice_bohr * atom_positions[:, i]
        @printf(io, "%d %12.6f %12.6f %12.6f %12.6f\n", n, charge, pos...)
    end

    for ix in 1:n_x
        for iy in 1:n_y
            for iz in 1:n_z
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
