using Printf: @printf
using Dates: now

# Cube format
# Specification from http://paulbourke.net/dataformats/cube/

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
