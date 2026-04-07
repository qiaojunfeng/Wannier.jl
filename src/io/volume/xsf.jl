# XSF format
# Specification from http://www.xcrysden.org/doc/XSF.html

export read_xsf, write_xsf

# TODO update
"""
    read_xsf(filename::AbstractString)

Read `xsf` file.

# Return
- `primvec`: `3 * 3`, Å, each column is a primitive lattice vector
- `convvec`: `3 * 3`, Å, each column is a conventional lattice vector
- `atoms`: `n_atoms` String, atomic symbols or numbers
- `atom_positions`: `3 * n_atoms`, Å, cartesian coordinates
- `rgrid`: [`RGrid`](@ref), grid on which `W` is defined
- `W`: `nx * ny * nz`, volumetric data

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

See also [`WannierIO.write_xsf`](@ref)
"""
function write_xsf(
        filename::AbstractString,
        lattice::AbstractMatrix,
        atom_positions::AbstractVector{V},
        atom_numbers::AbstractVector,
        rgrid::Union{RGrid, Nothing} = nothing,
        W::Union{AbstractArray{T, 3}, Nothing} = nothing,
    ) where {V <: AbstractVector, T <: Real}
    # Build a WannierIO.Xsf and delegate to the single-argument write_xsf API.
    # Convert atom numbers to string labels (WannierIO expects integer-parsable labels)
    atoms = string.(Int.(atom_numbers))

    # atom_positions are fractional coordinates w.r.t. `lattice`; convert to Cartesian
    atom_positions_cart = vec3.(frac_to_cart(lattice, atom_positions))

    if isnothing(rgrid) || isnothing(W)
        xsf = WannierIO.Xsf(mat3(lattice), mat3(lattice), atoms, atom_positions_cart, nothing, nothing, nothing, nothing, nothing, nothing)
        return WannierIO.write_xsf(filename, xsf)
    end

    O = origin(rgrid)
    spanvec = span_vectors(rgrid)
    # Build X, Y, Z fractional grids consistent with WannierIO expectations
    n_x, n_y, n_z = size(W)
    X = collect(range(0.0, 1.0, n_x))
    Y = collect(range(0.0, 1.0, n_y))
    Z = collect(range(0.0, 1.0, n_z))

    xsf = WannierIO.Xsf(
        mat3(lattice),
        mat3(lattice),
        atoms,
        atom_positions_cart,
        Vec3{Float64}(O),
        mat3(spanvec),
        X,
        Y,
        Z,
        W,
    )

    return WannierIO.write_xsf(filename, xsf)
end
