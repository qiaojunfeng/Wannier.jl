using Printf: @printf
using LinearAlgebra
import NearestNeighbors as NN

export get_bvectors

"""
    BVectorShells

Shells of bvectors.

The bvectors are sorted by norm such that equal-norm bvectors are grouped into one shell.

# Fields
- `recip_lattice`: `3 * 3`, each column is a reciprocal lattice vector
- `kpoints`: `3 * n_kpts`, in fractional coordinates
- `bvectors`: vectors of `3 * n_bvecs_per_shell`, in cartesian coordinates
- `weights`: vector of float, weights of each shell
- `multiplicities`: number of bvectors in each shell
- `n_bvecs`: total number of bvectors
"""
struct BVectorShells{T<:Real}
    # reciprocal lattice, 3 * 3, Å⁻¹ unit, each column is a lattice vector
    recip_lattice::Mat3{T}

    # kpoints array, fractional coordinates, 3 * n_kpts
    kpoints::Matrix{T}

    # bvectors of each shell, Cartesian! coordinates, Å⁻¹ unit
    # n_shells of 3 * multiplicity
    bvectors::Vector{Matrix{T}}

    # weight of each shell, length = n_shells
    weights::Vector{T}

    # multiplicity of each shell, length = n_shells
    multiplicities::Vector{Int}

    # number of shells
    n_shells::Int
end

"""
    BVectorShells(recip_lattice, kpoints, bvectors, weights)

Constructor of `BVectorShells`.

Only essential arguments are required, remaing fields of `BVectorShells`
are initialized accordingly.
This should be used instead of directly constructing `BVectorShells`.

# Arguments
- `recip_lattice`: `3 * 3`, each column is a reciprocal lattice vector
- `kpoints`: `3 * n_kpts`, in fractional coordinates
- `bvectors`: vectors of `3 * n_bvecs_per_shell`, in cartesian coordinates
- `weights`: vector of float, weights of each shell
"""
function BVectorShells(
    recip_lattice::Mat3{T},
    kpoints::Matrix{T},
    bvectors::Vector{Matrix{T}},
    weights::Vector{T},
) where {T<:Real}
    n_shells = length(bvectors)
    multiplicities = [size(bvectors[i], 2) for i in 1:n_shells]

    return BVectorShells{T}(
        recip_lattice, kpoints, bvectors, weights, multiplicities, n_shells
    )
end

function Base.show(io::IO, shells::BVectorShells)
    for i in 1:(shells.n_shells)
        @printf(io, "b-vector shell %3d    weight = %8.5f\n", i, shells.weights[i])
        vecs = shells.bvectors[i]
        for ib in axes(vecs, 2)
            ending = (i == shells.n_shells) && (ib == size(vecs, 2)) ? "" : "\n"
            @printf(io, "  %3d    %10.5f %10.5f %10.5f%s", ib, vecs[:, ib]..., ending)
        end
    end
end

"""
    BVectors

The bvectors for each kpoint.

# Fields
- `recip_lattice`: `3 * 3`, each column is a reciprocal lattice vector
- `kpoints`: `3 * n_kpts`, in fractional coordinates
- `bvectors`: `3 * n_bvecs`, in cartesian coordinates
- `weights`: `n_bvecs`, weights of each bvector
- `kpb_k`: k+b vectors at kpoint `k`, `k` -> `k + b`
    (index of periodically equivalent kpoint inside `recip_lattice`)
- `kpb_b`: `3 * n_bvecs * n_kpts`,
    displacements between k + b and its periodic image inside `recip_lattice`,
    such that k+b = `kpoints[:, kpb_k[ib, ik]] + kpb_b[:, ib, ik]` (in fractional)
- `n_kpts`: number of kpoints
- `n_bvecs`: total number of bvectors

!!! note

    In principle, we don't need to sort the bvectors for each kpoint, so that
    the bvectors have the same order as each kpoint.
    However, since `Wannier90` sort the bvectors, and the `mmn` file is written
    in that order, so we also sort in the same order as `Wannier90`.
"""
struct BVectors{T<:Real}
    # reciprocal lattice, 3 * 3, Å⁻¹ unit
    # each column is a reciprocal lattice vector
    recip_lattice::Mat3{T}

    # kpoints array, fractional coordinates, 3 * n_kpts
    kpoints::Matrix{T}

    # bvectors, Cartesian! coordinates, Å⁻¹ unit, 3 * n_bvecs
    bvectors::Matrix{T}

    # weight of each bvec, n_bvecs
    weights::Vector{T}

    # k+b vectors at kpoint k, n_bvecs * n_kpts
    # k -> k + b (index of periodically equivalent kpoint inside recip_lattice)
    kpb_k::Matrix{Int}

    # displacements between k + b and k + b wrapped around into the recip_lattice,
    # fractional coordinates, actually always integers since they are
    # the number of times to shift along each recip lattice.
    # 3 * n_bvecs * n_kpts, where 3 is [b_x, b_y, b_z]
    kpb_b::Array{Int,3}
end

function Base.getproperty(x::BVectors, sym::Symbol)
    if sym == :n_kpts
        return size(x.kpoints, 2)
    elseif sym == :n_bvecs
        return size(x.bvectors, 2)
    else
        # fallback to getfield
        getfield(x, sym)
    end
end

function Base.show(io::IO, bvectors::BVectors)
    println(io, "b-vectors:")
    @printf(io, "         [bx, by, bz] / Å⁻¹                weight\n")
    for i in 1:(bvectors.n_bvecs)
        v = bvectors.bvectors[:, i]
        w = bvectors.weights[i]
        ending = i < bvectors.n_bvecs ? "\n" : ""
        @printf(io, "%3d    %10.5f %10.5f %10.5f %10.5f%s", i, v..., w, ending)
    end
end

"""
    search_shells(kpoints, recip_lattice; atol=1e-6, max_shells=36)

Search bvector shells satisfing B1 condition.

# Arguments
- `kpoints`: fractional coordinates
- `recip_lattice`: each column is a reciprocal lattice vector

# Keyword Arguments
- `atol`: tolerance to select a shell (points having equal distances),
    equivalent to `Wannier90` input parameter `kmesh_tol`.
- `max_shells`: max number of nearest-neighbor shells,
    equivalent to `Wannier90` input parameter `search_shells`.
"""
function search_shells(
    kpoints::Matrix{T}, recip_lattice::Mat3{T}; atol::T=1e-6, max_shells::Int=36
) where {T<:Real}
    # Usually these "magic" numbers work well for normal recip_lattice.
    # Number of nearest-neighbors to be returned
    max_neighbors = 500
    # Max number of stencils in one shell
    max_multiplicity = 40
    # max_shells = round(Int, max_neighbors / max_multiplicity)

    # 1. Generate a supercell to search bvectors
    supercell, _ = make_supercell(kpoints)

    # To cartesian coordinates
    supercell_cart = similar(supercell)
    for ik in axes(supercell, 2)
        supercell_cart[:, ik] = recip_lattice * supercell[:, ik]
    end
    # use the 1st kpt to search bvectors, usually Gamma point
    kpt_orig = recip_lattice * kpoints[:, 1]

    # 2. KDTree to search nearest neighbors
    kdtree = NN.KDTree(supercell_cart)
    idxs, dists = NN.knn(kdtree, kpt_orig, max_neighbors, true)

    # activate debug info with: JULIA_DEBUG=Main julia
    # @debug "KDTree nearest neighbors" dists
    # @debug "KDTree nearest neighbors" idxs

    # 3. Arrange equal-distance kpoint indexes in layer of shells
    shells = Vector{Vector{Int}}()

    # The 1st result is the search point itself, dist = 0
    inb = 2
    ish = 1
    while (inb <= max_neighbors) && (ish <= max_shells)
        # use the 1st kpoint to find bvector shells & weights
        eqdist_idxs = findall(x -> isapprox(x, dists[inb]; atol=atol), dists)
        multi = length(eqdist_idxs)
        if multi >= max_multiplicity
            # skip large multiplicity shells
            inb += multi
            break
        end
        push!(shells, idxs[eqdist_idxs])
        inb += multi
        ish += 1
    end

    # 4. Get Cartesian coordinates vectors
    n_shells = length(shells)
    bvectors = Vector{Matrix{T}}(undef, n_shells)
    for ish in 1:n_shells
        kpb_cart = supercell_cart[:, shells[ish]]
        bvectors[ish] = kpb_cart .- kpt_orig
    end
    @debug "Found bvector shells" bvectors

    weights = zeros(T, 0)

    return BVectorShells(recip_lattice, kpoints, bvectors, weights)
end

"""
    are_parallel(A, B; atol=1e-6)

Check if the columns of matrix `A` and columns of matrix `B` are parallel.

# Arguments
- `A`: matrix
- `B`: matrix

# Keyword Arguments
- `atol`: tolerance to check parallelism
"""
function are_parallel(A::Matrix{T}, B::Matrix{T}; atol::T=1e-6) where {T<:Real}
    n_dim = size(A, 1)
    if n_dim != 3 || size(B, 1) != n_dim
        error("only support 3-vectors")
    end

    nc_A = size(A, 2)
    nc_B = size(B, 2)
    checkerboard = fill(false, nc_A, nc_B)

    for (i, c1) in enumerate(eachcol(A))
        for (j, c2) in enumerate(eachcol(B))
            p = cross(c1, c2)
            if all(isapprox.(0, p; atol=atol))
                checkerboard[i, j] = true
            end
        end
    end

    return checkerboard
end

"""
    delete_parallel(bvectors::Vector{Matrix{T}})

Remove shells having parallel bvectors.

# Arguments
- `bvectors`: vector of bvectors in each shell
"""
function delete_parallel(bvectors::Vector{Matrix{T}}) where {T<:Real}
    n_shells = length(bvectors)
    keep_shells = collect(1:n_shells)

    for ish in 2:n_shells
        for jsh in 1:(ish - 1)
            if !(jsh in keep_shells)
                continue
            end

            p = are_parallel(bvectors[jsh], bvectors[ish])
            if any(p)
                @debug "has parallel bvectors $jsh, $ish"
                filter!(s -> s != ish, keep_shells)
                break
            end
        end
    end
    new_bvectors = bvectors[keep_shells]

    @debug "keep shells" keep_shells
    @debug "After delete_parallel" [size(b, 2) for b in new_bvectors]' new_bvectors
    return new_bvectors
end

"""
    delete_parallel(shells::BVectorShells)

Remove shells having parallel bvectors.

# Arguments
- `shells`: `BVectorShells` containing bvectors in each shell
"""
function delete_parallel(shells::BVectorShells)
    bvectors = delete_parallel(shells.bvectors)
    return BVectorShells(shells.recip_lattice, shells.kpoints, bvectors, shells.weights)
end

"""
    compute_weights(bvectors::Vector{Matrix{T}}; atol=1e-6)

Try to guess bvector weights from MV1997 Eq. (B1).

The input bvectors are overcomplete vectors found during shell search, i.e. from `search_shells`.
This function tries to find the minimum number of bvector shells that
satisfy the B1 condition, and return the new `BVectorShells` and weights.

# Arguments
- `bvectors`: vector of bvectors in each shell

# Keyword Arguments
- `atol`: tolerance to satisfy B1 condition,
    equivalent to `Wannier90` input parameter `kmesh_tol`
"""
function compute_weights(bvectors::Vector{Matrix{T}}; atol::T=1e-6) where {T<:Real}
    n_shells = length(bvectors)
    n_shells == 0 && error("empty bvectors?")
    for i in 1:n_shells
        if size(bvectors[i], 1) != 3
            error("only support 3-vectors")
        end
    end

    # only compare the upper triangular part of bvec * bvec', 6 elements
    B = zeros(T, 6, n_shells)

    # return the upper triangular part of a matrix as a vector
    triu2vec(m) = m[triu!(trues(size(m)), 0)]
    # triu2vec(I) = [1 0 1 0 0 1]
    triu_I = triu2vec(diagm([1, 1, 1]))

    W = zeros(n_shells)

    # sigular value tolerance, to reproduce W90 behavior
    σ_atol = 1e-5

    keep_shells = zeros(Int, 0)
    ish = 1
    while ish <= n_shells
        push!(keep_shells, ish)
        B[:, ish] = triu2vec(bvectors[ish] * bvectors[ish]')
        # Solve equation B * W = triu_I
        # size(B) = (6, ishell), W is diagonal matrix of size ishell
        # B = U * S * V' -> W = V * S^-1 * U' * triu_I
        U, S, V = svd(B[:, keep_shells])
        @debug "S" ish S = S' keep_shells = keep_shells'
        if all(S .> σ_atol)
            W[keep_shells] = V * Diagonal(S)^-1 * U' * triu_I
            BW = B[:, keep_shells] * W[keep_shells]
            @debug "BW" ish BW = BW'
            if isapprox(BW, triu_I; atol=atol)
                break
            end
        else
            pop!(keep_shells)
        end
        ish += 1
    end
    if ish == n_shells + 1
        error("not enough shells to satisfy B1 condition")
    end
    n_shells = length(keep_shells)

    # resize
    weights = W[keep_shells]
    new_bvectors = bvectors[keep_shells]

    return new_bvectors, weights
end

"""
    compute_weights(shells::BVectorShells; atol=1e-6)

Try to guess bvector weights from MV1997 Eq. (B1).

# Arguments
- `shells`: `BVectorShells` containing bvectors in each shell

# Keyword Arguments
- `atol`: tolerance to satisfy B1 condition,
    equivalent to `Wannier90` input parameter `kmesh_tol`
"""
function compute_weights(shells::BVectorShells{T}; atol::T=1e-6) where {T<:Real}
    bvectors, weights = compute_weights(shells.bvectors; atol=atol)
    return BVectorShells(shells.recip_lattice, shells.kpoints, bvectors, weights)
end

"""
    check_b1(shells::BVectorShells; atol=1e-6)

Check completeness (B1 condition) of `BVectorShells`.

# Arguments
- `shells`: `BVectorShells` containing bvectors in each shell

# Keyword Arguments
- `atol`: tolerance, equivalent to `Wannier90` input parameter `kmesh_tol`
"""
function check_b1(shells::BVectorShells{T}; atol::T=1e-6) where {T<:Real}
    M = zeros(T, 3, 3)

    for ish in 1:(shells.n_shells)
        bvec = shells.bvectors[ish]
        M += shells.weights[ish] * bvec * bvec'
    end

    @debug "Bvector sum" M
    Δ = M - Matrix(I, 3, 3)
    # compare element-wise, to be consistent with W90
    if !all(isapprox.(Δ, 0; atol=atol))
        msg = "B1 condition is not satisfied\n"
        msg *= "  kmesh_tol = $atol\n"
        msg *= "  Δ = $(maximum(abs.(Δ)))\n"
        msg *= "  try increasing kmesh_tol?"
        error(msg)
    end

    println("Finite difference condition satisfied")
    println()
    return nothing
end

"""
    flatten_shells(shells::BVectorShells)

Flatten shell vectors into a matrix.

Return a tuple of `(bvecs, bvecs_weight)`, where
- `bvecs`: `3 * n_bvecs`
- `bvecs_weight`: `n_bvecs`
"""
function flatten_shells(shells::BVectorShells{T}) where {T<:Real}
    n_bvecs = sum(size(shells.bvectors[i], 2) for i in 1:(shells.n_shells))

    bvecs = zeros(T, 3, n_bvecs)
    bvecs_weight = zeros(T, n_bvecs)

    counter = 1
    for ish in 1:(shells.n_shells)
        multi = shells.multiplicities[ish]
        bvecs[:, counter:(counter + multi - 1)] = shells.bvectors[ish]
        bvecs_weight[counter:(counter + multi - 1)] .= shells.weights[ish]
        counter += multi
    end

    return bvecs, bvecs_weight
end

"""
    sort_supercell(translations, recip_lattice; atol=1e-8)

Sort supercell to fix the order of bvectors.

Both input and output `translations` are in fractional coordinates.

# Arguments
- `translations`: `3 * n_supercell` matrix, in fractional coordinates
- `recip_lattice`: each column is a reciprocal lattice vector

# Keyword Arguments
- `atol`: tolerance to compare bvectors,
    this is the same as what is hardcoded in `Wannier90`

!!! note

    This is used to reproduce `Wannier90` bvector order.
"""
function sort_supercell(
    translations::AbstractMatrix, recip_lattice::AbstractMatrix; atol::Real=1e-8
)
    n_cells = size(translations, 2)
    distances = zeros(eltype(recip_lattice), n_cells)

    for i in 1:n_cells
        distances[i] = norm(recip_lattice * translations[:, i])
    end

    # In W90, if the distances are degenerate, the distance which has larger index
    # is put at the front. Weird but I need to reproduce this.
    # So I reverse it, then do a stable sort in ascending order, so that this is a
    # "reversed" stable sort in ascending order.
    rev_distances = reverse(distances)
    # W90 atol = 1e-8, and need the equal sign
    lt(x, y) = x <= y - atol
    # The permutation is guaranteed to be stable
    perm = sortperm(rev_distances; lt=lt)

    idxs = collect(n_cells:-1:1)[perm]

    return translations[:, idxs]
end

"""
Find equivalent kpoint and displacement vector of bvectors `bvecs` at kpoint `k`.

all inputs in fractional coordinates.
"""
function _bvec_to_kb(bvecs::AbstractMatrix, k::AbstractVector, kpoints::AbstractMatrix)
    n_bvecs = size(bvecs, 2)

    kpts_equiv = zeros(Int, n_bvecs)
    b_equiv = zeros(Int, 3, n_bvecs)

    """Equivalent to periodic image?"""
    isequiv(v1, v2; atol=1e-6) = begin
        d = v1 - v2
        d -= round.(d)
        return all(isapprox.(d, 0; atol=atol))
    end

    for ib in 1:n_bvecs
        kpb = k + bvecs[:, ib]
        ik = findvector(isequiv, kpb, kpoints)
        kpts_equiv[ib] = ik
        b_equiv[:, ib] = round.(Int, kpb - kpoints[:, ik])
    end

    return kpts_equiv, b_equiv
end

"""
Sort bvectors specified by equivalent kpoint indices `k` and cell displacements `b`.

Sorting order:
1. length of bvectors: nearest k+b goes first, this is achieved by comparing
    the norm `bvecs_norm`.
2. supercell index: the supercell are already sorted by `sort_supercell`,
    which generates our input `translations`.
3. index of kpoint: the smaller index goes first, dictated by the input `kpoints`.

bvecs_norm: length of bvectors, cartesian norm.
k: index in `kpoints` for equivalent kpoint of bvectors
b: cell displacements of bvectors, fractional coordinates
translations: of supercell, fractional coordinates
"""
function _sort_kb(
    bvecs_norm::AbstractVector,
    k::AbstractVector{Int},
    b::AbstractMatrix{Int},
    translations::AbstractMatrix;
    atol::Real=1e-6,
)
    n_bvecs = length(k)

    # this is for comparing fractional coordinates, 1e-6 should be already safe
    isequiv(v1, v2) = isapprox(v1, v2; atol=1e-6)

    b_idx = zeros(Int, n_bvecs)
    for ib in 1:n_bvecs
        b_idx[ib] = findvector(isequiv, b[:, ib], translations)
    end

    lt(i, j) = begin
        if bvecs_norm[i] <= bvecs_norm[j] - atol
            return true
        elseif bvecs_norm[i] >= bvecs_norm[j] + atol
            return false
        else
            if b_idx[i] < b_idx[j]
                return true
            elseif b_idx[i] == b_idx[j]
                if k[i] < k[j]
                    return true
                else
                    return false
                end
            else
                return false
            end
        end
    end

    perm = sortperm(1:n_bvecs; lt=lt)

    return perm
end

"""
    sort_bvectors(shells::BVectorShells; atol=1e-6)

Sort bvectors in shells at each kpoints, to be consistent with `Wannier90`.

`Wannier90` use different order of bvectors at each kpoint,
in principle, this is not needed. However, the `mmn` file is
written in such order, so we need to sort bvectors and calculate
weights, since `nnkp` file has no section of weights.

# Arguments
- `shells`: `BVectorShells`

# Keyword Arguments
- `atol`: equivalent to `Wannier90` input parameter `kmesh_tol`
"""
function sort_bvectors(shells::BVectorShells{T}; atol::T=1e-6) where {T<:Real}
    kpoints = shells.kpoints
    recip_lattice = shells.recip_lattice
    n_kpts = size(kpoints, 2)

    # To sort bvectors for each kpoints, I need to
    # calculate distances of supercells to original cell.
    # I only need one kpoint at Gamma.
    _, translations = make_supercell(zeros(T, 3, 1))
    translations = sort_supercell(translations, recip_lattice)

    bvecs, bvecs_weight = flatten_shells(shells)
    n_bvecs = size(bvecs, 2)
    bvecs_frac = inv(recip_lattice) * bvecs
    bvecs_norm = [norm(bvecs[:, i]) for i in 1:n_bvecs]

    # find k+b indexes
    kpb_k = zeros(Int, n_bvecs, n_kpts)
    kpb_b = zeros(Int, 3, n_bvecs, n_kpts)
    # weight
    kpb_w = zeros(T, n_bvecs, n_kpts)

    for ik in 1:n_kpts
        k = kpoints[:, ik]
        # use fractional coordinates to compare
        k_equiv, b_equiv = _bvec_to_kb(bvecs_frac, k, kpoints)
        perm = _sort_kb(bvecs_norm, k_equiv, b_equiv, translations; atol=atol)

        kpb_k[:, ik] = k_equiv[perm]
        kpb_b[:, :, ik] = b_equiv[:, perm]
        kpb_w[:, ik] = bvecs_weight[perm]
    end

    # kpb_weight is redundant
    @assert sum(abs.(kpb_w .- bvecs_weight)) < 1e-6

    @debug "k+b k" kpb_k
    @debug "k+b b" kpb_b
    @debug "k+b weights" bvecs_weight

    return BVectors(recip_lattice, kpoints, bvecs, bvecs_weight, kpb_k, kpb_b)
end

"""
    get_bvectors(kpoints, recip_lattice; kmesh_tol=1e-6)

Generate and sort bvectors for all the kpoints.

# Arguments
- `kpoints`: `3 * n_kpts`, kpoints in fractional coordinates
- `recip_lattice`: `3 * 3`, columns are reciprocal lattice vectors

# Keyword Arguments
- `kmesh_tol`: equivalent to `Wannier90` input parameter `kmesh_tol`
"""
function get_bvectors(
    kpoints::Matrix{T}, recip_lattice::Mat3{T}; kmesh_tol::T=1e-6
) where {T<:Real}
    # find shells
    shells = search_shells(kpoints, recip_lattice; atol=kmesh_tol)
    shells = delete_parallel(shells)
    shells = compute_weights(shells; atol=kmesh_tol)
    show(shells)
    println("\n")
    check_b1(shells; atol=kmesh_tol)

    # generate bvectors for each kpoint
    bvectors = sort_bvectors(shells; atol=kmesh_tol)

    return bvectors
end

"""
    index_bvector(kpb_k, kpb_b, k1, k2, b)

Given bvector `b` connecting kpoints `k1` and `k2`, return the index of the bvector `ib`.

This is a reverse search of bvector index if you only know the two kpoints `k1` and `k2`, and
the connecting displacement vector `b`.

# Arguments
- `kpb_k`: `n_bvecs * n_kpts`, k+b kpoints at `k1`
- `kpb_b`: `3 * n_bvecs * n_kpts`, displacement vector for k+b bvectors at `k1`
- `k1`: integer, index of kpoint `k1`
- `k2`: integer, index of kpoint `k2`
- `b`: vector of 3 integer, displacement vector from `k1` to `k2`
"""
function index_bvector(
    kpb_k::Matrix{Int}, kpb_b::Array{Int,3}, k1::Int, k2::Int, b::AbstractVector{Int}
)
    n_bvecs = size(kpb_k, 1)

    for ib in 1:n_bvecs
        if kpb_k[ib, k1] == k2 && kpb_b[:, ib, k1] == b
            return ib
        end
    end

    return error("No neighbors found, k1 = $(k1), k2 = $(k2), b = $(b)")
end

function index_bvector(bvectors::BVectors, k1, k2, b)
    return index_bvector(bvectors.kpb_k, bvectors.kpb_b, k1, k2, b)
end

"""
    get_bvectors_nearest(kpoints, recip_lattice; kmesh_tol=1e-6)

Generate and sort bvectors for all the kpoints.

# Arguments
- `kpoints`: `3 * n_kpts`, kpoints in fractional coordinates
- `recip_lattice`: `3 * 3`, columns are reciprocal lattice vectors

# Keyword Arguments
- `kmesh_tol`: equivalent to `Wannier90` input parameter `kmesh_tol`
"""
function get_bvectors_nearest(kpoints::Matrix{T}, recip_lattice::Mat3{T}) where {T<:Real}
    n_kx, n_ky, n_kz = get_kgrid(kpoints)
    δx, δy, δz = 1 / n_kx, 1 / n_ky, 1 / n_kz

    # only 6 nearest neighbors
    n_bvecs = 6
    bvecs_frac = zeros(T, 3, n_bvecs)
    bvecs_frac[:, 1] = [δx, 0, 0]
    bvecs_frac[:, 2] = [-δx, 0, 0]
    bvecs_frac[:, 3] = [0, δy, 0]
    bvecs_frac[:, 4] = [0, -δy, 0]
    bvecs_frac[:, 5] = [0, 0, δz]
    bvecs_frac[:, 6] = [0, 0, -δz]

    # just a fake weight
    bvecs_weight = ones(T, n_bvecs)

    # generate bvectors for each kpoint
    n_kpts = size(kpoints, 2)
    kpb_k = zeros(Int, n_bvecs, n_kpts)
    kpb_b = zeros(Int, 3, n_bvecs, n_kpts)

    for ik in 1:n_kpts
        k = @view kpoints[:, ik]
        # use fractional coordinates to compare
        k_equiv, b_equiv = _bvec_to_kb(bvecs_frac, k, kpoints)
        kpb_k[:, ik] = k_equiv
        kpb_b[:, :, ik] = b_equiv
    end

    bvecs = recip_lattice * bvecs_frac
    return BVectors(recip_lattice, kpoints, bvecs, bvecs_weight, kpb_k, kpb_b)
end
