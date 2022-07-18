using Printf: @printf
using LinearAlgebra
import NearestNeighbors as NN

@doc raw"""
Shells of bvectors.
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

    # multiplicity of each shell, length = num_shells
    multiplicities::Vector{Int}

    # number of shells
    n_shells::Int
end

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

function pprint(shells::BVectorShells)
    for i in 1:(shells.n_shells)
        @printf("b-vector shell %3d    weight = %8.5f\n", i, shells.weights[i])
        vecs = shells.bvectors[i]
        for ib in 1:size(vecs, 2)
            @printf("  %3d    %10.5f %10.5f %10.5f\n", ib, vecs[:, ib]...)
        end
    end
    # print a blank line to separate the following stdout
    println()
    return nothing
end

@doc raw"""
The bvectors for each kpoint, sorted in the same order as the W90 nnkp file.
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

    # k+b vectors, k -> k + b (index of equivalent kpt in the recip_cell), n_bvecs * n_kpts
    kpb_k::Matrix{Int}

    # displacements between k + b and k + b wrapped around into the recip_cell,
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

function pprint(bvectors::BVectors)
    println("b-vectors:")
    @printf("         [bx, by, bz] / Å⁻¹                weight\n")
    for i in 1:(bvectors.n_bvecs)
        v = bvectors.bvectors[:, i]
        w = bvectors.weights[i]
        @printf("%3d    %10.5f %10.5f %10.5f %10.5f\n", i, v..., w)
    end
    # print a blank line to separate the following stdout
    println()
    return nothing
end

@doc raw"""
Make a supercell of kpoints by translating it along 3 directions.
Input and returned kpoints are in fractional coordinates.

repeat: number of repetitions along ±x, ±y, ±z directions, on output
there are (2*repeat + 1)^3 cells.
"""
function make_supercell(kpoints::Matrix{T}, replica::Int=5) where {T<:Real}
    size(kpoints, 1) ≉ 3 && error("kpoints must be 3 * n_kpts")
    n_kpts = size(kpoints, 2)

    n_cell = (2 * replica + 1)^3

    supercell = Matrix{T}(undef, 3, n_cell * n_kpts)
    translations = Matrix{Int}(undef, 3, n_cell * n_kpts)

    counter = 1
    for ix in (-replica):replica
        for iy in (-replica):replica
            for iz in (-replica):replica
                for ik in 1:n_kpts
                    supercell[:, counter] = kpoints[:, ik] + [ix, iy, iz]
                    translations[:, counter] = [ix, iy, iz]
                    counter += 1
                end
            end
        end
    end

    return supercell, translations
end

@doc raw"""
Search bvector shells satisfing B1 condition.

kpoints: fractional coordinates
recip_lattice: each column is a reciprocal lattice vector
atol: tolerance to select a shell (points having equal distances),
    equivalent to W90 `kmesh_tol`.
max_shells: max number of nearest-neighbor shells,
    equivalent to W90 `search_shells`.
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
    for ik in 1:size(supercell, 2)
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

@doc raw"""
Check if the columns of matrix A and columns of matrix B are parallel.
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

@doc raw"""
Remove shells having parallel bvectors
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

function delete_parallel(shells::BVectorShells)
    bvectors = delete_parallel(shells.bvectors)
    return BVectorShells(shells.recip_lattice, shells.kpoints, bvectors, shells.weights)
end

@doc raw"""
Try to guess bvector weights from MV1997 Eq. (B1).

input bvectors are overcomplete vectors found during shell search.
atol: tolerance to satisfy B1 condition, equivalent to W90 `kmesh_tol`
return: bvectors and weights satisfing B1 condition.
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

function compute_weights(shells::BVectorShells{T}; atol::T=1e-6) where {T<:Real}
    bvectors, weights = compute_weights(shells.bvectors; atol=atol)
    return BVectorShells(shells.recip_lattice, shells.kpoints, bvectors, weights)
end

@doc raw"""
Check completeness (B1 condition) of bvectors.

atol: tolerance, equivalent to W90 `kmesh_tol`
"""
function check_b1(shells::BVectorShells{T}; atol::T=1e-6) where {T<:Real}
    M = zeros(T, 3, 3)

    for ish in 1:(shells.n_shells)
        bvec = shells.bvectors[ish]
        M += shells.weights[ish] * bvec * bvec'
    end

    @debug "Bvector sum" M
    if !isapprox(M, I; atol=atol)
        error("B1 condition is not satisfied")
    end

    println("Finite difference condition satisfied")
    println()
    return nothing
end

@doc raw"""
flatten shell vectors into a matrix
return:
bvecs: 3 * num_bvecs
bvecs_weight: num_bvecs
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
Sort supercell to fix the order of bvectors.
To reproduce W90 bvec order.

translations: 3 x n_supercell matrix of fractional coordinates
output translations: fractional coordinates
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

"""Find vector in the columns of a matrix"""
function findvector(predicate::Function, v::AbstractVector, M::AbstractMatrix)
    for (i, col) in enumerate(eachcol(M))
        predicate(v, col) && return i
    end
    error("$v not found in array!")
    return nothing
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

@doc raw"""
Sort bvectors in shells for each kpoints, to be consistent with wannier90

Wannier90 use different order of bvectors for each kpoint,
in principle, this is not needed. However, the Mmn file is
written in such order, so I need to sort bvectors and calculate
weights, since nnkp file has no section of weights.
atol: equivalent to W90 `kmesh_tol`
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

@doc raw"""
Generate bvectors for all the kpoints.
"""
function get_bvectors(
    kpoints::Matrix{T}, recip_lattice::Mat3{T}; kmesh_tol::T=1e-6
) where {T<:Real}
    # find shells
    shells = search_shells(kpoints, recip_lattice; atol=kmesh_tol)
    shells = delete_parallel(shells)
    shells = compute_weights(shells; atol=kmesh_tol)
    pprint(shells)
    check_b1(shells; atol=kmesh_tol)

    # generate bvectors for each kpoint
    bvectors = sort_bvectors(shells; atol=kmesh_tol)

    return bvectors
end
