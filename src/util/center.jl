import NearestNeighbors as NN

"""
    find_nearests(point, search_neighbors, lattice, atom_positions)

Find nearest-atom (including its periodic images) to a `point`.

Usually `point` is the WF center, so the function returns nearest atom to WF center.

# Arguments
- `point`: vector of `3` floats, fractional coordinates w.r.t. lattice
- `search_neighbors`: number of nearest-neighbors to be returned
- `lattice`: `3 * 3`, each column is a lattice vector
- `atom_positions`: `3 * n_atoms`, each column is fractional coordinate ∈ `[0, 1)`
"""
function find_nearests(
    point::AbstractVector{T},
    search_neighbors::R,
    lattice::AbstractMatrix{T},
    atom_positions::Vector,
) where {T<:Real,R<:Integer}
    @assert length(point) == 3
    @assert search_neighbors > 0
    @assert size(lattice) == (3, 3)

    point_cart = lattice * point
    periodic_atoms, periodic_translations = make_supercell(atom_positions)
    periodic_atoms_cart = map(p -> lattice * p, periodic_atoms)

    kdtree = NN.KDTree(periodic_atoms_cart)
    idxs, distances = NN.knn(kdtree, point_cart, search_neighbors, true)

    # find corresponding atoms in the original unit cell
    translations = periodic_translations[idxs]
    atom_wrapped_cart = periodic_atoms_cart[idxs] .- map(t -> lattice * t, translations)

    # back to fractional coordinates so that the comparison is independent of lattice length
    atom_wrapped = map(a -> inv(lattice) * a, atom_wrapped_cart)
    # find indexes of atom in atom_positions
    indexes = map(x -> findfirst(isapprox(x; atol=1e-6), atom_positions), atom_wrapped)

    return distances, indexes, translations
end

"""
    find_nearest_atom(centers, lattice, atom_positions)

Find nearest atom for each WF center.

# Arguments
- `centers`: `3 * n_wann`, in fractional coordinates
- `lattice`: `3 * 3`, each column is a lattice vector
- `atom_positions`: `3 * n_atoms`, each column is fractional coordinate ∈ `[0, 1)`

# Example
```julia
wout = read_wout("silicon.wout")
points = inv(wout.lattice) * wout.centers  # to fractional
find_nearest_atom(points, wout.lattice, wout.atom_positions)
```
"""
function find_nearest_atom(
    centers::Vector, lattice::AbstractMatrix{T}, atom_positions::Vector{Vec3{T}}
) where {T<:Real}
    n_wann = length(centers)
    distances = zeros(T, n_wann)
    indexes = zeros(Int, n_wann)
    translations = zeros(Vec3{Int}, n_wann)

    for iw in 1:n_wann
        d, i, t = find_nearests(centers[iw], 1, lattice, atom_positions)
        distances[iw] = d[1]
        indexes[iw] = i[1]
        translations[iw] = t[1]
    end
    return distances, indexes, translations
end

"""
    wrap_centers(centers, lattice)

Wrap around centers back to unit cell at origin.

# Arguments
- `centers`:: `3 * n_wann`, in Cartesian coordiantes
- `lattice`:: `3 * 3`, each column is a lattice vector
"""
function wrap_centers(centers::Vector, lattice::AbstractMatrix)
    inv_lattice = inv(lattice)
    return map(centers) do c
        # mod: back to home cell
        lattice * mod.(inv_lattice * c, 1)
    end
end
