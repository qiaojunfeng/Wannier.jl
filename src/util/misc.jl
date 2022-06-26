import LinearAlgebra as LA
import NearestNeighbors as NN


function get_projectability(A)
    num_bands, num_wann, num_kpts = size(A)
    proj = zeros(num_bands, num_kpts)
    for ik = 1:num_kpts
        p = A[:, :, ik] * A[:, :, ik]'
        proj[:, ik] = real(LA.diag(p))
    end
    return proj
end


@doc raw"""
Find nearest-atom to a point (usually it is the center of a Wannier function)
reduced_coord: atoms, point are relative to unit_cell, in reduced coordinates

unit_cell: 3 * 3, each column is a lattice vector
atoms: 3 * num_atoms, each column is a coordinate
point: 3
search_neighbors: number of nearest-neighbors to be returned
"""
function find_nearests(unit_cell::T, atoms::T, point::Vector{Float64}
    ; reduced_coord::Bool=true, search_neighbors::Int=5) where {T<:Union{Matrix{Float64},LA.Adjoint{Float64,Matrix{Float64}}}}

    @assert size(unit_cell) == (3, 3)
    @assert size(atoms, 1) == 3
    @assert length(point) == 3
    @assert search_neighbors > 0

    inv_cell = inv(unit_cell)
    if reduced_coord
        atoms_cart = unit_cell * atoms
        point_cart = unit_cell * point
        atoms_reduced = atoms
        point_reduced = point
    else
        atoms_cart = atoms
        point_cart = point
        atoms_reduced = inv_cell * atoms
        point_reduced = inv_cell * point
    end

    reapeated_atoms_reduced, supercell_idx = generate_supercell(atoms_reduced)
    repeated_atoms_cart = unit_cell * reapeated_atoms_reduced

    kdtree = NN.KDTree(repeated_atoms_cart)
    idxs, dists = NN.knn(kdtree, point_cart, search_neighbors, true)

    # find corresponding atoms in the original unit cell
    translations = supercell_idx[:, idxs]
    atom_unit_cell = repeated_atoms_cart[:, idxs] - unit_cell * translations
    function findvec(mat, v, dims)
        iseqv(v1, v2) = all(isapprox.(vec(v1), vec(v2); atol=1e-6))
        found = mapslices(elem -> iseqv(elem, v), mat; dims=dims)
        return findfirst(dropdims(found; dims=1))
    end
    idxs_unit_cell = [findvec(atoms_cart, atom_unit_cell[:, i], 1) for i = 1:size(atom_unit_cell, 2)]

    return dists, idxs_unit_cell, translations
end
