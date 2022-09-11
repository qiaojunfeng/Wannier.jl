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
    cube = WannierIO.read_cube(filename)
    rgrid = RGrid(cube.span_vectors, cube.origin, cube.X, cube.Y, cube.Z)
    return (; cube.atom_positions, cube.atom_numbers, rgrid, cube.W)
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

    return WannierIO.write_cube(filename, all_pos, all_atoms, O, spanvec, W)
end
