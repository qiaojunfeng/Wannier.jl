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
    xsf = WannierIO.read_xsf(filename)
    rgrid = nothing

    if !isnothing(xsf.W)
        rgrid = RGrid(xsf.span_vectors, xsf.origin, xsf.X, xsf.Y, xsf.Z)
    end

    return (; xsf.primvec, xsf.convvec, xsf.atoms, xsf.atom_positions, rgrid, xsf.W)
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
    return WannierIO.write_xsf(
        filename, lattice, atom_positions, atom_numbers, O, spanvec, W
    )
end
