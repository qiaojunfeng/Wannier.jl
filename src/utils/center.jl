using NearestNeighbors: KDTree, knn

export find_neighbors, find_nearest_neighbor

"""
    $(SIGNATURES)

Find several neighboring atoms (including its periodic images) to several `points`.

Usually `points` are WF centers, so the function returns from the 1st to the
n-th nearest neighboring atoms to WF centers.

# Arguments
- `points`: a vector of length `n_points`, each element is a 3-vector for
    fractional coordinates w.r.t. lattice
- `n_neighbors`: number of nearest-neighbors to be returned
- `lattice`: `3 * 3`, each column is a lattice vector
- `atom_positions`: length-`n_atoms` vector, each element is a 3-vector for
    fractional coordinate ∈ `[0, 1)`

# Return
- `distances`: length-`n_points` vector, each element is a length-`n_neighbors`
    vector of distances to the nearest neighbors
- `indices`: length-`n_points` vector, each element is a length-`n_neighbors`
    vector of indexes of the nearest neighbors in `atom_positions`
- `translations`: length-`n_points` vector, each element is a length-`n_neighbors`
    vector of translations from the atoms in `atom_positions` to the real
    periodically repeated atoms which are the actual nearest neighbors of `points`
"""
function find_neighbors(
    points::AbstractVector{T},
    n_neighbors::Integer,
    lattice::AbstractMatrix,
    atom_positions::AbstractVector,
) where {T<:AbstractVector}
    @assert n_neighbors > 0
    @assert size(lattice) == (3, 3)

    periodic_atoms, periodic_translations = make_supercell(atom_positions)
    # force to Vec3 for `KDTree`
    periodic_atoms_cart = map(p -> Vec3(lattice * p), periodic_atoms)
    kdtree = KDTree(periodic_atoms_cart)
    points_cart = map(p -> Vec3(lattice * p), points)
    idxs, distances = knn(kdtree, points_cart, n_neighbors, true)

    # find corresponding atoms in the original unit cell
    translations = [periodic_translations[idx] for idx in idxs]
    # find translations from original atom_positions to the found periodic atoms
    indices = map(enumerate(idxs)) do (i, idx)
        # in fractional coordinates so that the comparison is independent of lattice length
        atom_wrapped = periodic_atoms[idx] .- translations[i]
        # find indexes of atom in atom_positions
        map(v -> findfirst(isapprox(v; atol=1e-6), atom_positions), atom_wrapped)
    end

    return distances, indices, translations
end

"""
    $(SIGNATURES)

Find several neighboring atoms (including its periodic images) to a single `point`.

See also [`find_neighbors`](@ref) for mutliple input points.
"""
function find_neighbors(
    point::AbstractVector{T},
    n_neighbors::Integer,
    lattice::AbstractMatrix,
    atom_positions::AbstractVector,
) where {T<:Real}
    @assert length(point) == 3
    distances, indices, translations = find_neighbors(
        [point], n_neighbors, lattice, atom_positions
    )
    # unwrap the result
    return distances[1], indices[1], translations[1]
end

"""
    $(SIGNATURES)

Find the nearest atom for several WF centers.

# Arguments
- `centers`: a vector, each element is a 3-vector for fractional coordinates
- `lattice`: `3 * 3`, each column is a lattice vector
- `atom_positions`: a vector, each element is a fractional coordinate ∈ `[0, 1)`

# Example
```julia
wout = read_wout("silicon.wout")
centers = inv(wout.lattice) * wout.centers  # to fractional
find_nearest_neighbor(centers, wout.lattice, wout.atom_positions)
```
"""
function find_nearest_neighbor(
    centers::AbstractVector{T}, lattice::AbstractMatrix, atom_positions::AbstractVector
) where {T<:AbstractVector}
    @assert length(centers) > 0 "Empty centers"
    # unwrap the result
    distances, indices, translations = find_neighbors(centers, 1, lattice, atom_positions)
    distances = map(d -> d[1], distances)
    indices = map(i -> i[1], indices)
    translations = map(t -> t[1], translations)
    return distances, indices, translations
end

"""
    $(SIGNATURES)

Wrap centers back to unit cell at origin.

# Arguments
- `centers`:: a vector, each element is a 3-vector for Cartesian coordiantes
- `lattice`:: `3 * 3`, each column is a lattice vector

# Return
- centers inside unit cell, in Cartesian coordinates

!!! note

    Here the input `centers` is in Cartesian coordinates; for fractional
    coordinates, you just need to remove the the integer part of the coordinates.
"""
function wrap_centers(centers::AbstractVector, lattice::AbstractMatrix)
    inv_lattice = inv(lattice)
    return map(centers) do c
        # mod: back to home cell
        lattice * mod.(inv_lattice * c, 1)
    end
end
