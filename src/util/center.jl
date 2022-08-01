import NearestNeighbors as NN

@doc raw"""
Find nearest-atom (including its periodic images) to a point
(usually it is the center of a Wannier function).

point: 3, fractional coordinates relative to lattice.
search_neighbors: number of nearest-neighbors to be returned
lattice: 3 * 3, each column is a lattice vector
atom_positions: 3 * num_atoms, each column is fractional coordinate ∈ [0, 1)
"""
function find_nearests(
    point::AbstractVector{T},
    search_neighbors::R,
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
) where {T<:Real,R<:Integer}
    @assert length(point) == 3
    @assert search_neighbors > 0
    @assert size(lattice) == (3, 3)
    @assert size(atom_positions, 1) == 3

    point_cart = lattice * point

    periodic_atoms, periodic_translations = make_supercell(atom_positions)
    periodic_atoms_cart = lattice * periodic_atoms

    kdtree = NN.KDTree(periodic_atoms_cart)
    idxs, distances = NN.knn(kdtree, point_cart, search_neighbors, true)

    # find corresponding atoms in the original unit cell
    translations = periodic_translations[:, idxs]
    atom_wrapped_cart = periodic_atoms_cart[:, idxs] - lattice * translations

    # back to fractional coordinates so that the comparison is independent of lattice length
    atom_wrapped = inv(lattice) * atom_wrapped_cart
    # find indexes of atom in atom_positions
    iseqv(v1, v2) = all(isapprox.(vec(v1), vec(v2); atol=1e-6))
    f(v) = findvector(iseqv, v, atom_positions)
    indexes = f.([v for v in eachcol(atom_wrapped)])

    return distances, indexes, translations
end

"""
Find nearest atom of WFs.

centers: 3 * n_wann, in fractional coordinates
lattice: 3 * 3, each column is a lattice vector
atom_positions: 3 * n_atoms, each column is fractional coordinate ∈ [0, 1)

Example usage:
```
wout = read_wout("silicon.wout")
points = inv(wout.lattice) * wout.centers  # to fractional
find_nearest_atom(points, wout.lattice, wout.atom_positions)
```
"""
function find_nearest_atom(
    centers::AbstractMatrix{T},
    lattice::AbstractMatrix{T},
    atom_positions::AbstractMatrix{T},
) where {T<:Real}
    n_wann = size(centers, 2)
    distances = zeros(T, n_wann)
    indexes = zeros(Int, n_wann)
    translations = zeros(Int, 3, n_wann)

    for iw in 1:n_wann
        d, i, t = find_nearests(centers[:, iw], 1, lattice, atom_positions)
        distances[iw] = d[1]
        indexes[iw] = i[1]
        translations[:, iw] = t[:, 1]
    end
    return distances, indexes, translations
end

"""
Wrap around centers back to unit cell at origin

centers:: 3 * n_wann, in Cartesian! coordiantes
lattice:: 3 * 3, each column is a lattice vector
"""
function wrap_centers(
    centers::AbstractMatrix{T}, lattice::AbstractMatrix{T}
) where {T<:Real}
    centers_frac = inv(lattice) * centers
    # back to home cell
    centers_frac = mod.(centers_frac, 1)
    return lattice * centers_frac
end
