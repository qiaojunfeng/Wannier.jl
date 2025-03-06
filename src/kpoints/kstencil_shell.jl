using NearestNeighbors: knn, KDTree

"""
    $(TYPEDEF)

Shells of b-vectors.

The ``\\mathbf{b}``-vectors are are vectors connecting a kpoint to its neighboring
kpoints. To find a stencil for approximating finite differences, the neighboring
kpoints are sorted by their distance to the original kpoint, such that equal-distance
``\\mathbf{b}``-vectors are grouped into one shell.

# Fields
$(FIELDS)
"""
struct KspaceStencilShells{T<:Real}
    """reciprocal lattice vectors, 3 * 3, each column is a reciprocal lattice
    vector in Å⁻¹ unit"""
    recip_lattice::Mat3{T}

    """number of kpoints along three reciprocal lattice vectors"""
    kgrid_size::Vec3{Int}

    """kpoint fractional coordinates, length-`n_kpoints` vector of `Vec3`.
    Should be a uniformlly-spaced kpoint grid."""
    kpoints::Vector{Vec3{T}}

    """dengeneracy (i.e. the number of bvectors) of each shell,
    length-`n_shells` vector of integers"""
    n_degens::Vector{Int}

    """bvectors of each shell, Cartesian coordinates in Å⁻¹ unit,
    length-`n_shells` vector, each element is a length-`n_degens[i_shell]`
    vector of `Vec3`.

    Here we use Cartesian coordinates because it is easier to compute the
    completeness condition (MV1997 Eq. (B1)).
    """
    bvectors::Vector{Vector{Vec3{T}}}

    """bvector weight of each shell, length-`n_shells` vector, Å² unit."""
    bweights::Vector{T}
end

n_kpoints(shells::KspaceStencilShells) = length(shells.kpoints)
"""number of b-vector shells"""
n_shells(shells::KspaceStencilShells) = length(shells.n_degens)
n_bvectors(shells::KspaceStencilShells) = sum(shells.n_degens)
reciprocal_lattice(shells::KspaceStencilShells) = shell.recip_lattice

"""
    $(SIGNATURES)

Convenience constructor of `KspaceStencilShells`, auto set `n_degens`.

# Arguments
See the fields of [`KspaceStencilShells`](@ref) struct.
"""
function KspaceStencilShells(recip_lattice, kgrid_size, kpoints, bvectors, bweights)
    n_degens = [length(bvecs) for bvecs in bvectors]
    return KspaceStencilShells(
        recip_lattice, Vec3(kgrid_size), kpoints, n_degens, bvectors, bweights
    )
end

function KspaceStencilShells(
    recip_lattice, kgrid_size, kpoints;
    atol=default_w90_kmesh_tol(),
)
    # find shells
    shells = search_shells(recip_lattice, kgrid_size, kpoints; atol)
    keep_shells = check_parallel(shells)
    shells = delete_shells(shells, keep_shells)

    keep_shells, bweights = compute_bweights(shells; atol)
    shells = delete_shells(shells, keep_shells)
    shells.bweights .= bweights

    # Γ-point calculation only keep half of the bvectors
    if all(kgrid_size .== 1)
        shells = delete_shells_Γ(shells)
    end

    check_completeness(shells; atol)
    return shells
end

function Base.show(io::IO, ::MIME"text/plain", shells::KspaceStencilShells)
    nshells = n_shells(shells)
    @printf(io, "                 [bx, by, bz] (Å⁻¹)\n")
    for (ish, (bvecs, w)) in enumerate(zip(shells.bvectors, shells.bweights))
        n = isempty(bvecs) ? NaN : norm(bvecs[1])
        @printf(
            io, "b-vector shell %3d:    norm = %8.5f (Å⁻¹)   weight = %8.5f (Å²)", ish, n, w
        )
        for (ib, bvec) in enumerate(shells.bvectors[ish])
            @printf(io, "\n%3d    %11.6f %11.6f %11.6f", ib, bvec...)
        end
        ish != nshells && println(io)
    end
end

"""
    $(SIGNATURES)

Search bvector shells satisfing completeness condition.

# Arguments
- `recip_lattice`: each column is a reciprocal lattice vector
- `kpoints`: fractional coordinates

# Keyword Arguments
- `atol`: tolerance to select a shell (points having equal distances)
- `max_shells`: max number of nearest-neighbor shells

# Return
- a `KspaceStencilShells` struct, note the `bweights` are not computed yet, all zeros!

!!! note
    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`
    - `max_shells` should be set to wannier90's input parameter `search_shells`
"""
function search_shells(
    recip_lattice::Mat3,
    kgrid_size::AbstractVector,
    kpoints::AbstractVector;
    atol=default_w90_kmesh_tol(),
    max_shells=default_w90_bvectors_search_shells(),
)
    # Usually these "magic" numbers work well for normal recip_lattice.
    # Number of nearest-neighbors to be returned
    max_neighbors = 500
    # Max number of stencils in one shell
    max_degens = 40
    # max_shells = round(Int, max_neighbors / max_degens)

    # 1. Generate a supercell to search bvectors
    supercell, _ = make_supercell(kpoints)
    # To cartesian coordinates
    supercell_cart = map(supercell) do cell
        recip_lattice * cell
    end
    # use the 1st kpt to search bvectors, usually Γ point
    kpt_orig = recip_lattice * kpoints[1]

    # 2. KDTree to search nearest neighbors
    kdtree = KDTree(supercell_cart)
    idxs, dists = knn(kdtree, kpt_orig, max_neighbors, true)
    # activate debug info with: JULIA_DEBUG=Main julia
    # @debug "KDTree nearest neighbors" dists
    # @debug "KDTree nearest neighbors" idxs

    # 3. Arrange equal-distance kpoint indices in layer of shells
    shells = Vector{Vector{Int}}()

    # The 1st result is the search point itself, dist = 0
    inb = 2  # index of neighbors
    ish = 1  # index of shells
    while (inb <= max_neighbors) && (ish <= max_shells)
        # use the 1st kpoint to find bvector shells & bweights
        eqdist_idxs = findall(x -> isapprox(x, dists[inb]; atol), dists)
        degen = length(eqdist_idxs)
        if degen >= max_degens
            # skip large-degeneracy shells
            inb += degen
            break
        end
        push!(shells, idxs[eqdist_idxs])
        inb += degen
        ish += 1
    end

    # 4. Get Cartesian-coordinate vectors
    T = eltype(recip_lattice)
    bvectors = map(shells) do idxs
        kpb_cart = supercell_cart[idxs]
        return map(v -> Vec3{T}(v .- kpt_orig), kpb_cart)
    end
    @debug "Found bvector shells" bvectors

    bweights = zeros(T, length(shells))
    return KspaceStencilShells(recip_lattice, kgrid_size, kpoints, bvectors, bweights)
end

"""
    $(SIGNATURES)

Check if the columns of matrix `A` and columns of matrix `B` are parallel.

# Arguments
- `A`: vector, each element is a vector, often `Vec3`
- `B`: similar to `A`

# Keyword Arguments
- `atol`: tolerance to check parallelism.

# Return
- `checkerboard`: boolean matrix, `checkerboard[i, j]` is `true`
    if `A[i]` and `B[j]` are parallel

!!! note
    Wannier90 uses a constant `1e-6` as `atol`, thus here its default is set to
    `1e-6` as well to reproduce the same result.
"""
function are_parallel(
    A::AbstractVector, B::AbstractVector; atol=default_w90_bvectors_check_parallel_atol()
)
    n_A = length(A)
    n_B = length(B)

    checkerboard = fill(false, n_A, n_B)

    for (i, a) in enumerate(A)
        for (j, b) in enumerate(B)
            p = cross(a, b)
            if all(isapprox.(0, p; atol))
                checkerboard[i, j] = true
            end
        end
    end

    return checkerboard
end

"""
    $(SIGNATURES)

Check if shells having parallel bvectors.

# Arguments
- `bvectors`: vector of bvectors in each shell

# Keyword Arguments
- `atol`: tolerance to check parallelism, see [`are_parallel`](@ref)

# Return
- `keep_shells`: indices of shells that do not have parallel bvectors

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's internal constant `1e-6`.
"""
function check_parallel(
    bvectors::Vector{Vector{Vec3{T}}}; atol=default_w90_bvectors_check_parallel_atol()
) where {T}
    nshells = length(bvectors)
    keep_shells = collect(1:nshells)

    for ish in 2:nshells
        for jsh in 1:(ish - 1)
            if !(jsh in keep_shells)
                continue
            end

            p = are_parallel(bvectors[jsh], bvectors[ish]; atol)
            if any(p)
                @debug "has parallel bvectors between shells $jsh $ish"
                filter!(s -> s != ish, keep_shells)
                break
            end
        end
    end
    return keep_shells
end

"""
    $(SIGNATURES)

Check if shells having parallel bvectors.

# Arguments
- `shells`: `KspaceStencilShells` containing bvectors in each shell
"""
function check_parallel(shells::KspaceStencilShells)
    return check_parallel(shells.bvectors)
end

function delete_shells(bvectors::Vector{Vector{Vec3{T}}}, keep_shells) where {T}
    return bvectors[keep_shells]
end

"""
    $(SIGNATURES)

Remove shells.

# Arguments
- `keep_shells`: indices of shells to keep
"""
function delete_shells(shells::KspaceStencilShells, keep_shells)
    bvectors = delete_shells(shells.bvectors, keep_shells)
    bweights = shells.bweights[keep_shells]
    return KspaceStencilShells(
        shells.recip_lattice, shells.kgrid_size, shells.kpoints, bvectors, bweights
    )
end

"""
    $(SIGNATURES)

Delete negetive bvectors for Γ-point calculation.

Since bvectors are symmetric, this removes half of the bvectors.
"""
function delete_shells_Γ(shells::KspaceStencilShells)
    bvectors = map(shells.bvectors) do bvecs  # for each shell
        bvecs_new = filter(v -> all(v .>= 0), bvecs)
        if length(bvecs_new) != length(bvecs)//2
            error("Non-symmetric bvectors for Γ-point calculation: ", bvecs)
        end
        bvecs_new
    end
    bweights = [2w for w in shells.bweights]
    return KspaceStencilShells(
        shells.recip_lattice, shells.kgrid_size, shells.kpoints, bvectors, bweights
    )
end

"""
    $(SIGNATURES)

Try to guess bvector bweights from MV1997 Eq. (B1).

The input bvectors are overcomplete vectors found during shell search, i.e., from
[`search_shells`](@ref). This function tries to find the minimum number of bvector
shells that satisfy the B1 condition, and return the new `KspaceStencilShells` and bweights.

# Arguments
- `bvectors`: vector of bvectors in each shell

# Keyword Arguments
- `atol`: tolerance to satisfy B1 condition

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`
"""
function compute_bweights(
    bvectors::Vector{Vector{Vec3{T}}}; atol=default_w90_kmesh_tol()
) where {T}
    nshells = length(bvectors)
    @assert nshells > 0 "empty bvectors"

    # only compare the upper triangular part of bvec * bvec', 6 elements
    B = zeros(T, 6, nshells)

    # return the upper triangular part of a matrix column-by-column as a vector
    # e.g., triu2vec(I) = [1 0 1 0 0 1]
    triu2vec(m::AbstractMatrix) = m[triu!(trues(size(m)), 0)]
    triu_I = triu2vec(diagm([1, 1, 1]))

    # weight of each shell
    W = zeros(nshells)

    # sigular value tolerance, to reproduce W90 behavior
    σ_atol = default_w90_bvectors_singular_value_atol()

    keep_shells = zeros(Int, 0)
    ish = 1
    while ish <= nshells
        push!(keep_shells, ish)
        # to 3 * n_degens[ish] matrix
        b = reduce(hcat, bvectors[ish])
        B[:, ish] = triu2vec(b * b')
        # Solve equation B * W = triu_I
        # size(B) = (6, n_shells), W is diagonal matrix of size n_shells
        # B = U * S * V' -> W = V * S^-1 * U' * triu_I
        U, S, V = svd(B[:, keep_shells])
        @debug "S" ish S = S' keep_shells = keep_shells'
        if all(S .> σ_atol)
            W[keep_shells] = V * inv(Diagonal(S)) * U' * triu_I
            BW = B[:, keep_shells] * W[keep_shells]
            @debug "BW" ish BW = BW'
            if isapprox(BW, triu_I; atol)
                break
            end
        else
            pop!(keep_shells)
        end
        ish += 1
    end
    if ish == nshells + 1
        error("not enough shells to satisfy B1 condition")
    end

    bweights = W[keep_shells]
    return keep_shells, bweights
end

"""
    $(SIGNATURES)

Try to guess bvector bweights from MV1997 Eq. (B1).

# Arguments
- `shells`: `KspaceStencilShells` containing bvectors in each shell

# Keyword Arguments
- `atol`: tolerance to satisfy B1 condition

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`
"""
function compute_bweights(shells::KspaceStencilShells; atol=default_w90_kmesh_tol())
    return compute_bweights(shells.bvectors; atol)
end

"""
    $(SIGNATURES)

Check completeness (B1 condition) of `KspaceStencilShells`.

# Arguments
- `shells`: `KspaceStencilShells` containing bvectors in each shell

# Keyword Arguments
- `atol`: floating point tolerance

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`
"""
function check_completeness(
    shells::KspaceStencilShells{T}; atol=default_w90_kmesh_tol()
) where {T}
    M = zeros(T, 3, 3)

    for (bvecs, w) in zip(shells.bvectors, shells.bweights)
        bvec = reduce(hcat, bvecs)
        M += w * bvec * bvec'
    end

    @debug "Bvector sum" M
    Δ = M - Matrix(I, 3, 3)
    # compare element-wise, to be consistent with W90
    if !all(isapprox.(Δ, 0; atol))
        error("""b-vector completeness condition not satisfied
                 atol = $atol
                 Δ = $(maximum(abs.(Δ)))
                 try increasing atol?""")
    end

    @info "b-vector completeness condition satisfied"
    return nothing
end

"""
    $(SIGNATURES)

Unwrap nested shell vectors into a flattened vector.

# Return
- `bvectors`: length-`n_bvectors` vector, each element is a `Vec3`
- `bweights`: length-`n_bvectors` vector of bweights for each bvector
"""
function flatten_shells(shells::KspaceStencilShells{T}) where {T}
    nbvecs = n_bvectors(shells)

    bvectors = zeros(Vec3{T}, nbvecs)
    bweights = zeros(T, nbvecs)

    counter = 1
    for (bvecs, w, degen) in zip(shells.bvectors, shells.bweights, shells.n_degens)
        bvectors[counter:(counter + degen - 1)] = bvecs
        bweights[counter:(counter + degen - 1)] .= w
        counter += degen
    end

    return bvectors, bweights
end
