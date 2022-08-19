using Printf: @printf
using Dates: now

# XSF format
# Specification from http://www.xcrysden.org/doc/XSF.html

"""
Read xsf file.

All outputs in cartesian coordinates, bohr unit.
"""
function read_xsf(filename::AbstractString)
    io = open(filename)

    primvec = nothing
    convvec = nothing
    atoms = nothing
    positions = nothing
    origin = nothing
    spanvec = nothing
    W = nothing

    while !eof(io)
        line = strip(readline(io))
        if isempty(line) || startswith(line, '#')
            continue
        end

        if occursin("CRYSTAL", line)
            # primitive lattice, each column is a lattice vector
            @assert strip(readline(io)) == "PRIMVEC"
            primvec = zeros(Float64, 3, 3)
            for i in 1:3
                line = strip(readline(io))
                primvec[:, i] = parse.(Float64, split(line))
            end
            # conventional lattice, each column is a lattice vector
            @assert strip(readline(io)) == "CONVVEC"
            convvec = zeros(Float64, 3, 3)
            for i in 1:3
                line = strip(readline(io))
                convvec[:, i] = parse.(Float64, split(line))
            end
            # read atom positions
            @assert strip(readline(io)) == "PRIMCOORD"
            line = strip(readline(io))
            n_atom = parse(Int, split(line)[1])
            atoms = Vector{String}(undef, n_atom)
            # each column is a position vector
            positions = zeros(Float64, 3, n_atom)
            for i in 1:n_atom
                line = strip(readline(io))
                # might be element label, or atomic number
                atoms[i] = split(line)[1]
                positions[:, i] = parse.(Float64, split(line)[2:4])
            end
        elseif occursin("BEGIN_BLOCK_DATAGRID_3D", line)
            comment = strip(readline(io))
            line = strip(readline(io))
            # I only read the 1st data grid, others are ignored
            @assert startswith(line, "BEGIN_DATAGRID_3D")
            identifier = chopprefix(line, "BEGIN_DATAGRID_3D_")
            ngx, ngy, ngz = parse.(Int, split(strip(readline(io))))
            origin = parse.(Float64, split(strip(readline(io))))
            # spanning vectors
            spanvec = zeros(Float64, 3, 3)
            for i in 1:3
                line = strip(readline(io))
                spanvec[:, i] = parse.(Float64, split(line))
            end
            # column-major
            W = zeros(Float64, ngx, ngy, ngz)
            idx = 1
            while idx <= ngx * ngy * ngz
                line = split(strip(readline(io)))
                ncol = length(line)
                W[idx:(idx + ncol - 1)] = parse.(Float64, line)
                idx += ncol
            end
            @assert occursin("END_DATAGRID_3D", strip(readline(io)))
            @assert strip(readline(io)) == "END_BLOCK_DATAGRID_3D"
        end
    end

    return (
        primvec=primvec,
        convvec=convvec,
        atoms=atoms,
        positions=positions,
        origin=origin,
        spanvec=spanvec,
        W=W,
    )
end

"""
Write xsf file.

atom_positions: 3 x n_atoms, fractional coordinates
atom_numbers: n_atoms, atomic numbers
lattice: 3 x 3, Å, each column is a lattice vector
n_voxels: 3, number of voxels along three lattice vectors
X: fractional coordinates of W along a1
Y: fractional coordinates of W along a2
Z: fractional coordinates of W along a3
W: nx x ny x nz, volumetric data
"""
function write_xsf(
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

    @info "Writing xsf file: " filename
    io = open(filename, "w")

    # header
    @printf(io, "# Created by Wannier.jl %s\n", string(now()))

    @printf(io, "CRYSTAL\n")
    @printf(io, "PRIMVEC\n")
    @printf(io, "%12.7f %12.7f %12.7f\n", lattice[:, 1]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", lattice[:, 2]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", lattice[:, 3]...)
    @printf(io, "CONVVEC\n")
    @printf(io, "%12.7f %12.7f %12.7f\n", lattice[:, 1]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", lattice[:, 2]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", lattice[:, 3]...)
    @printf(io, "PRIMCOORD\n")
    @printf(io, "%d 1\n", n_atoms)
    for i in 1:n_atoms
        pos = lattice * atom_positions[:, i]
        @printf(io, "%d %12.7f %12.7f %12.7f\n", atom_numbers[i], pos...)
    end

    @printf(io, "\n")
    @printf(io, "BEGIN_BLOCK_DATAGRID_3D\n")
    @printf(io, "3D_field\n")
    @printf(io, "BEGIN_DATAGRID_3D_UNKNOWN\n")

    ngx, ngy, ngz = size(W)
    @printf(io, "%d %d %d\n", ngx, ngy, ngz)
    # origin
    origin = lattice * [X[1], Y[1], Z[1]]
    @printf(io, "%12.7f %12.7f %12.7f\n", origin...)
    ax1 = lattice * [X[end] - X[1], 0, 0]
    ax2 = lattice * [0, Y[end] - Y[1], 0]
    ax3 = lattice * [0, 0, Z[end] - Z[1]]
    @printf(io, "%12.7f %12.7f %12.7f\n", ax1...)
    @printf(io, "%12.7f %12.7f %12.7f\n", ax2...)
    @printf(io, "%12.7f %12.7f %12.7f\n", ax3...)

    # column-major
    ncol = 0
    for k in 1:ngz
        for j in 1:ngy
            for i in 1:ngx
                @printf(io, " %13.5e", W[i, j, k])
                ncol += 1
                if ncol == 6
                    @printf(io, "\n")
                    ncol = 0
                end
            end
        end
    end
    ncol != 0 && @printf(io, "\n")
    @printf(io, "END_DATAGRID_3D\n")
    @printf(io, "END_BLOCK_DATAGRID_3D\n")

    close(io)
    return nothing
end
