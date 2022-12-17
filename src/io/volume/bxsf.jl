export read_bxsf, write_bxsf

"""
    read_bxsf(filename::AbstractString)

Read `bxsf` file.

# Return
- `rgrid`: [`RGrid`](@ref), grid on which `E` is defined
- `fermi_energy`: in eV
- `E`: `n_bands * n_kx * n_ky * n_kz`, energy eigenvalues
"""
function read_bxsf(filename::AbstractString)
    bxsf = WannierIO.read_bxsf(filename)
    rgrid = nothing

    if !isnothing(bxsf.E)
        rgrid = RGrid(bxsf.span_vectors, bxsf.origin, bxsf.X, bxsf.Y, bxsf.Z)
    end

    return (; rgrid, bxsf.fermi_energy, bxsf.E)
end

"""
    write_bxsf(filename, lattice, atom_positions, atom_numbers, rgrid, W)

Write `bxsf` file.

# Arguments
- `rgrid`: `RGrid`
- `fermi_energy`: in eV
- `E`: `n_bands * n_kx * n_ky * n_kz`, energy eigenvalues

This is a more user-friendly version. The `rgrid` contains the information of the
grid origin and spanning vectors.

See also [`WannierIO.write_bxsf`](@ref)
"""
function write_bxsf(
    filename::AbstractString, rgrid::RGrid, fermi_energy::T, E::AbstractArray{T,4}
) where {T<:Real}
    O = origin(rgrid)
    spanvec = span_vectors(rgrid)
    return WannierIO.write_bxsf(filename, fermi_energy, O, spanvec, E)
end
