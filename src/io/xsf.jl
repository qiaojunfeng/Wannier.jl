using Printf: @printf
using Dates: now

# XSF format
# Specification from http://www.xcrysden.org/doc/XSF.html

export read_xsf, write_xsf

"""
    read_xsf(filename::AbstractString)

Read `xsf` file.

All outputs in cartesian coordinates, Å unit.

!!! note

    Only support reading 1 datagrid in `BLOCK_DATAGRID_3D`.
"""
function read_xsf(filename::AbstractString)
    io = open(filename)

    primvec = nothing
    convvec = nothing
    atoms = nothing
    atom_positions = nothing
    origin = nothing
    span_vectors = nothing
    rgrid = nothing
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
            atom_positions = zeros(Float64, 3, n_atom)
            for i in 1:n_atom
                line = strip(readline(io))
                # might be element label, or atomic number
                atoms[i] = split(line)[1]
                atom_positions[:, i] = parse.(Float64, split(line)[2:4])
            end
        elseif occursin("BEGIN_BLOCK_DATAGRID_3D", line)
            comment = strip(readline(io))
            line = strip(readline(io))
            # I only read the 1st data grid, others are ignored
            @assert startswith(line, "BEGIN_DATAGRID_3D")
            # identifier = chopprefix(line, "BEGIN_DATAGRID_3D_")
            ngx, ngy, ngz = parse.(Int, split(strip(readline(io))))
            origin = parse.(Float64, split(strip(readline(io))))
            # spanning vectors
            span_vectors = zeros(Float64, 3, 3)
            for i in 1:3
                line = strip(readline(io))
                span_vectors[:, i] = parse.(Float64, split(line))
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

    if !isnothing(W)
        n_x, n_y, n_z = size(W)
        # fractional w.r.t. span_vectors
        O = inv(span_vectors) * origin
        X = range(0, 1, n_x) .+ O[1]
        Y = range(0, 1, n_y) .+ O[2]
        Z = range(0, 1, n_z) .+ O[3]
        Xg, Yg, Zg = ndgrid(X, Y, Z)
        rgrid = RGrid(span_vectors, Xg, Yg, Zg)
    end

    return (; primvec, convvec, atoms, atom_positions, rgrid, W)
end

"""
    write_xsf(filename, lattice, atom_positions, atom_numbers, origin, span_vectors, W)

Write `xsf` file.

# Arguments
- `lattice`: `3 * 3`, Å, each column is a lattice vector
- `atom_positions`: `3 * n_atoms`, fractional coordinates
- `atom_numbers`: `n_atoms`, atomic numbers
- `origin`: `3`, Å, origin of the grid
- `span_vectors`: `3 * 3`, Å, each column is a spanning vector
- `W`: `nx * ny * nz`, volumetric data

See also [`write_xsf(filename, lattice, atom_positions, atom_numbers, rgrid, W)`]
(@ref write_xsf(filename, lattice, atom_positions, atom_numbers, rgrid, W))
"""
function write_xsf(
    filename::AbstractString,
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{Int},
    origin::AbstractVector{T},
    span_vectors::AbstractMatrix{T},
    W::AbstractArray{T,3},
) where {T<:Real}
    n_atoms = length(atom_numbers)
    size(atom_positions, 2) == n_atoms || error("incompatible n_atoms")
    size(lattice) == (3, 3) || error("incompatible lattice")
    size(span_vectors) == (3, 3) || error("incompatible span_vectors")

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

    n_x, n_y, n_z = size(W)
    @printf(io, "%d %d %d\n", n_x, n_y, n_z)
    @printf(io, "%12.7f %12.7f %12.7f\n", origin...)
    @printf(io, "%12.7f %12.7f %12.7f\n", span_vectors[:, 1]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", span_vectors[:, 2]...)
    @printf(io, "%12.7f %12.7f %12.7f\n", span_vectors[:, 3]...)

    # column-major
    ncol = 0
    for k in 1:n_z
        for j in 1:n_y
            for i in 1:n_x
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

"""
    write_xsf(filename, lattice, atom_positions, atom_numbers, rgrid, W)

Write `xsf` file.

# Arguments
- `lattice`: `3 * 3`, Å, each column is a lattice vector
- `atom_positions`: `3 * n_atoms`, fractional coordinates
- `atom_numbers`: `n_atoms`, atomic numbers
- `rgrid`: `RGrid`
- `W`: `nx * ny * nz`, volumetric data

This is a more user-friendly version. The `rgrid` contains the information of the
grid origin and spanning vectors.

See also [`write_xsf(filename, lattice, atom_positions, atom_numbers, origin, span_vectors, W)`]
(@ref write_xsf(filename, lattice, atom_positions, atom_numbers, origin, span_vectors, W))
"""
function write_xsf(
    filename::AbstractString,
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
    atom_numbers::AbstractVector{Int},
    rgrid::RGrid,
    W::AbstractArray{T,3},
) where {T<:Real}
    O = origin(rgrid)
    spanvec = span_vectors(rgrid)
    return write_xsf(filename, lattice, atom_positions, atom_numbers, O, spanvec, W)
end
