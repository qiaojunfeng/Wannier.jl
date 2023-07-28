export generate_kspace_stencil, index_bvector

"""
    $(TYPEDEF)

The kpoint grid stencil (i.e., ``\\mathbf{b}``-vectors) at each kpoint.

# Fields
$(FIELDS)

!!! note

    In principle, we don't need to sort the bvectors for each kpoint, so that
    the bvectors have the same ordering at each kpoint.
    However, since wannier90 sorts the bvectors, if we want to read the `mmn`
    file, we need to make sure we sort the bvectors in the same ordering as
    wannier90.
"""
struct KspaceStencil{T<:Real}
    """reciprocal lattice vectors, 3 * 3, each column is a reciprocal lattice
    vector in Å⁻¹ unit"""
    recip_lattice::Mat3{T}

    """number of kpoints along three reciprocal lattice vectors"""
    kgrid_size::Vec3{Int}

    """kpoint fractional coordinates, length-`n_kpoints` vector of `Vec3`.
    Should be a uniformlly-spaced kpoint grid."""
    kpoints::Vector{Vec3{T}}

    """\\mathbf{b}-vectors, length-`n_bvectors` vector, each element is a `Vec3`
    in Cartesian coordinates, in Å⁻¹ unit"""
    bvectors::Vector{Vec3{T}}

    """weight of each bvector, length-`n_bvectors` vector, Å² unit"""
    bweights::Vector{T}

    """Indices of kpoints that are periodic images of \\mathbf{k+b} vectors.
    the real ``\\mathbf{k+b}`` vector has periodic image ``\\mathbf{k^\\prime}``
    kpoint inside the reciprocal lattice, the `kpb_k[ik][ib]` is the index
    of this ``\\mathbf{k^\\prime}`` kpoint of the `kpoints` variable, i.e.,
    the true ``\\mathbf{b}``-vector can be retrieved by
    `b = kpoints[kpb_k[ik][ib]] + kpb_G[ik][ib] - kpoints[ik]`.
    length-`n_kpoints` vector, each element is a length-`n_bvectors` vector
    of integers. `kpb` is the abbreviation for `k plus b`."""
    kpb_k::Vector{Vector{Int}}

    """displacement vector between the true ``\\mathbf{k+b}`` vector and its
    periodic image inside reciprocal lattice.
    The `k + b` is wrapped around into the reciprocal lattice, and by adding
    this `kpb_G` we restores the true `k+b` vector which can be outside of
    reciprocal lattice. See also the definition of `kpb_k`.

    In fractional coordinates, actually always integers since they are
    multiples of reciprocal lattice vectors.
    Length-`n_kpoints` vector, each element is a length-`n_bvectors` vector,
    then each element is a `Vec3` for the x, y, z components of the
    shifting vector."""
    kpb_G::Vector{Vector{Vec3{Int}}}
end

function KspaceStencil(recip_lattice, kgrid_size, kpoints, bvectors, bweights, kpb_k, kpb_G)
    return KspaceStencil(
        Mat3(recip_lattice), Vec3(kgrid_size), kpoints, bvectors, bweights, kpb_k, kpb_G
    )
end

n_kpoints(kstencil::KspaceStencil) = length(kstencil.kpoints)
n_bvectors(kstencil::KspaceStencil) = length(kstencil.bvectors)
reciprocal_lattice(kstencil::KspaceStencil) = kstencil.recip_lattice

function Base.show(io::IO, ::MIME"text/plain", kstencil::KspaceStencil)
    show_recip_lattice(io, kstencil.recip_lattice)
    println(io)
    @printf(io, "b-vectors:       [bx, by, bz] (Å⁻¹)          norm (Å⁻¹)  weight (Å²)")
    for (i, (v, w)) in enumerate(zip(kstencil.bvectors, kstencil.bweights))
        @printf(io, "\n%3d    %11.6f %11.6f %11.6f %11.6f %11.6f", i, v..., norm(v), w)
    end
    println(io, "\n")
    @printf(io, "kgrid_size  =  %d %d %d\n", kstencil.kgrid_size...)
    @printf(io, "n_kpoints   =  %d\n", n_kpoints(kstencil))
    @printf(io, "n_bvectors  =  %d", n_bvectors(kstencil))
end

"""
    $(SIGNATURES)

Compare two `KspaceStencil` objects.
"""
function Base.isapprox(a::KspaceStencil, b::KspaceStencil; kwargs...)
    for f in propertynames(a)
        va = getfield(a, f)
        vb = getfield(b, f)

        if va isa Vector
            all(isapprox.(va, vb; kwargs...)) || return false
        else
            isapprox(va, vb; kwargs...) || return false
        end
    end
    return true
end

"""
    $(SIGNATURES)

Compute bvector bweights from MV1997 Eq. (B1).

# Arguments
- `bvectors`: a **flattened** vector of bvectors, assuming the bvectors are
    already correct but only need to compute the bweights. e.g., after parsing
    the `mmn` or `nnkp` file, the bvector themselves are already known, but
    the bweights are still missing.

# Keyword Arguments
- `atol`: tolerance to satisfy B1 condition

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`
"""
function compute_bweights(bvectors::Vector{Vec3{T}}; atol=default_w90_kmesh_tol()) where {T}
    # assume the bvectors are correct: they should be able to be nested into
    # a shell structure
    bvectors_norm = map(norm, bvectors)
    # the input bvectors can be unordered, so I sort them to be safe, then
    # fold them into shells which are ordered by the norm of bvectors
    perm = sortperm(bvectors_norm)

    # permuted bvectors so that they are ordered by norm
    bvectors_sorted = bvectors[perm]
    # nest bvectors into equal-norm shells
    bvectors_nested = Vector{Vector{Vec3{T}}}(undef, 0)
    # indices of bvectors_sort in each shell
    shells = Vector{Vector{Int}}(undef, 0)

    nbvecs = length(bvectors)
    ib = 1
    bvectors_norm_sorted = bvectors_norm[perm]
    while ib <= nbvecs
        # find bvectors having equal norm
        b_norm = bvectors_norm_sorted[ib]
        eqnorm_idxs = findall(x -> isapprox(x, b_norm; atol), bvectors_norm_sorted)
        push!(bvectors_nested, bvectors_sorted[eqnorm_idxs])
        push!(shells, eqnorm_idxs)
        ib += length(eqnorm_idxs)
    end
    keep_shells, bweights = compute_bweights(bvectors_nested; atol)
    @assert keep_shells == 1:length(bvectors_nested) "compute_bweights should be " *
        " idempotent to the bvectors, maybe the input bvectors are not complete?"

    # now remap bweights to the original bvector order
    # mappings: index of bvectors_sorted -> index of shell
    mappings = [fill(i, length(sh)) for (i, sh) in enumerate(shells)]
    # flatten mappings
    mappings = vcat(mappings...)
    return map(perm) do p
        bweights[mappings[p]]
    end
end

"""
    $(SIGNATURES)

Sort supercell to fix the order of bvectors.

Both input and output `translations` are in fractional coordinates.

# Arguments
- `translations`: length-`n_supercell` vector, each element is a 3-vector for
    translations in fractional coordinates
- `recip_lattice`: each column is a reciprocal lattice vector in Å⁻¹ unit

# Keyword Arguments
- `atol`: tolerance to compare bvectors

# Return
- `translations`: sorted supercell translations

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's internal constant `1e-8`.
"""
function sort_supercell(
    translations::AbstractVector,
    recip_lattice::AbstractMatrix;
    atol=default_w90_bvectors_sort_supercell_atol(),
)
    n_cells = length(translations)
    distances = map(translations) do t
        norm(recip_lattice * t)
    end

    # In W90, if the distances are degenerate, the distance which has larger index
    # is put at the front. Weird but I need to reproduce this.
    # So I reverse it, then do a stable sort in ascending order, so that this is a
    # "reversed" stable sort in ascending order.
    rev_distances = reverse(distances)
    # W90 atol = 1e-8, and need the equal sign
    lt(x, y) = x <= y - atol
    # The permutation is guaranteed to be stable
    perm = sortperm(rev_distances; lt)

    idxs = collect(n_cells:-1:1)[perm]

    # The returned translations are in ascending order of distances w.r.t. to
    # origin. For degenerate distances, the larger index (in the original
    # translations) is put at the front.
    return translations[idxs]
end

"""
    $(SIGNATURES)

Find equivalent kpoints and displacement vectors of bvectors.

# Arguments
- `bvectors_frac`: length-`n_bvectors` vector, each element is a b-vectors in
    fractional coordinates
- `k`: fractional coodinates of a single kpoint, which can be outside of [0, 1);
    therefore, its `G ≠ [0, 0, 0]`.
- `kpoints`: length-`n_kpoints` vector of fractional coordinates of all kpoints

# Return
- `kpb_k`: length-`n_bvectors` vector of (shifted) `k+b` kpoint indices
- `kpb_G`: length-`n_bvectors` vector of displacement vectors to reach `k+b`
    kpoints, in fractional coordinates

!!! note

    All inputs in fractional coordinates for better floating-point comparisons.
"""
function bvectors_to_kpb(bvectors_frac::AbstractVector, k::Vec3, kpoints::AbstractVector)
    nbvecs = length(bvectors_frac)

    kpb_k = zeros(Int, nbvecs)
    kpb_G = zeros(Vec3{Int}, nbvecs)

    """
    Check if two vectors are equivalent apart from some integers.

    W90 internally use `kmesh_tol` to compare the norm of (v1 - v2) which are in
    Cartesian coordinates; while here we compare directly the fractional
    coordinates so it's safer (does not depend on the actual length of
    recip_lattice), thus I choose a (somewhat arbitrary) constant `1e-6`.
    """
    isequiv(v1, v2; atol=1e-6) = begin
        d = v1 - v2
        d -= round.(d)
        return all(isapprox.(d, 0; atol))
    end

    for ib in 1:nbvecs
        kpb = k + bvectors_frac[ib]
        ik = findfirst(k -> isequiv(kpb, k), kpoints)
        kpb_k[ib] = ik
        kpb_G[ib] = round.(Int, kpb - kpoints[ik])
    end

    return kpb_k, kpb_G
end

"""
    $(SIGNATURES)

Sort bvectors (at a single kpoint) specified by equivalent kpoint indices `iks`
and cell displacements `Gs`.

Sorting order:
1. length of bvectors: nearest k+b goes first, this is achieved by comparing
    the bvector norms `b_norms`.
2. supercell index: the supercell are already sorted by [`sort_supercell`](@ref),
    which generates our input `translations`.
3. index of kpoint: the smaller index goes first, dictated by the input `kpoints`.

# Arguments
- `b_norms`: norm of bvectors, Cartesian in Å⁻¹ unit.
- `iks`: indices in `kpoints` for equivalent kpoints of `k+b` vectors
- `Gs`: cell displacements to reach `k+b` vectors, fractional coordinates
- `translations`: of supercell, fractional coordinates

# Return
- `perm`: permutation of `1:n_bvectors` such that the bvectors at this kpoint
    are sorted in the order described above.
    e.g., if on input the bvectors are stored in a `bs::Vector{Vec3}`,
    then `bs[perm]` is the sorted bvectors.

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`
"""
function sort_kpb(
    b_norms::AbstractVector,
    iks::AbstractVector{Int},
    Gs::AbstractVector{Vec3{Int}},
    translations::AbstractVector{Vec3{Int}};
    atol=default_w90_kmesh_tol(),
)
    # this is for comparing fractional coordinates, 1e-6 should be already safe
    isequiv(v1, v2) = isapprox(v1, v2; atol=1e-6)

    G_idxs = map(Gs) do G
        findfirst(t -> isequiv(t, G), translations)
    end

    lt(i, j) = begin
        if b_norms[i] <= b_norms[j] - atol
            return true
        elseif b_norms[i] >= b_norms[j] + atol
            return false
        else
            if G_idxs[i] < G_idxs[j]
                return true
            elseif G_idxs[i] == G_idxs[j]
                if iks[i] < iks[j]
                    return true
                else
                    return false
                end
            else
                return false
            end
        end
    end

    perm = sortperm(1:length(iks); lt)
    # @debug "sort_kpb" perm
    return perm
end

"""
    $(SIGNATURES)

Sort bvectors in shells at each kpoint, to reproduce the same ordering as
wannier90.

Wannier90 uses different ordering of bvectors at each kpoint. In principle,
this is not needed. However, the `mmn` file is written in such order, so we need
to sort bvectors and calculate bweights, since `nnkp` file has no section of bweights.

# Arguments
- `shells`: `KspaceStencilShells`

# Keyword Arguments
- `atol`: floating point tolerance

!!! note

    To reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`
"""
function sort_bvectors(
    shells::KspaceStencilShells{T}; atol=default_w90_kmesh_tol()
) where {T}
    kpoints = shells.kpoints
    recip_lattice = shells.recip_lattice

    # To sort bvectors for each kpoints, I need to calculate distances of
    # supercells to original cell. I only need one kpoint at Γ.
    _, translations = make_supercell([Vec3{T}(0, 0, 0)])
    # Becareful, the `atol` of `sort_supercell` is different from `kmesh_tol`!
    # If using `kmesh_tol`, the dataset"graphene" nnkp file test will fail.
    translations = sort_supercell(translations, recip_lattice)

    bvectors, bweights = flatten_shells(shells)
    nbvecs = length(bvectors)
    bvectors_norm = norm.(bvectors)
    inv_recip_lattice = inv(recip_lattice)
    bvectors_frac = map(bvectors) do b
        inv_recip_lattice * b
    end

    # find k+b indices
    nkpts = length(kpoints)
    kpb_k = [zeros(Int, nbvecs) for _ in 1:nkpts]
    kpb_G = [zeros(Vec3{Int}, nbvecs) for _ in 1:nkpts]

    for (ik, kpt) in enumerate(kpoints)
        # use fractional coordinates for comparisons
        ik_equiv, G_equiv = bvectors_to_kpb(bvectors_frac, kpt, kpoints)
        perm = sort_kpb(bvectors_norm, ik_equiv, G_equiv, translations; atol)

        kpb_k[ik] = ik_equiv[perm]
        kpb_G[ik] = G_equiv[perm]
        # since small-length bvectors go first, their bweights should not change
        @assert iszero(bweights - bweights[perm]) "bvector bweights should not change"
    end

    @debug "k+b k" kpb_k
    @debug "k+b G" kpb_G
    @debug "k+b bweights" bweights

    # Reset the order of `bvectors` by using the ordering of the 1st kpoint.
    # In principle, this is not needed -- our `kpb_k` and `kpb_G` already have
    # the same ordering as wannier90.
    # However, this additional step can make sure the `KspaceStencil.bvectors` are
    # exactly the same as the wannier90 wout file, so we can directly test
    # against the wout file.
    bvectors = map(zip(kpb_k[1], kpb_G[1])) do (ikpb, G)
        recip_lattice * (kpoints[ikpb] + G - kpoints[1])
    end

    return KspaceStencil(
        recip_lattice, shells.kgrid_size, kpoints, bvectors, bweights, kpb_k, kpb_G
    )
end

"""Abstract type for b-vector generation algorithms"""
abstract type KspaceStencilAlgorithm end

"""Generate b-vectors for first-order finite difference, i.e., MV1997 Eq. (B1)"""
struct FirstOrderKspaceStencil <: KspaceStencilAlgorithm end

"""Generate b-vectors for 6 (cubic, ±x, ±y, ±z) nearest neighbors"""
struct CubicNearestKspaceStencil <: KspaceStencilAlgorithm end

"""
    generate_kspace_stencil()

Generate bvectors for all the kpoints.

# Arguments
- `kpoints`: kpoints in fractional coordinates

# Keyword Arguments
- `atol`: floating point tolerance

# Return
- a [`KspaceStencil`](@ref) struct

!!! tip

    The default algorithm is `FirstOrderKspaceStencil`, which generates bvectors
    same as wannier90, exactly in the same ordering.
    To fully reproduce wannier90's behavior,
    - `atol` should be set to wannier90's input parameter `kmesh_tol`

    The `CubicNearestKspaceStencil` generates bvectors containing only 6 nearest
    neighbors (cubic case). This is useful for [`parallel_transport`](@ref),
    since we only need overlap matrices between 6 nearest neighbors, and some
    times the `FirstOrderKspaceStencil` does not contain those 6 nearest neighbors.
"""
function generate_kspace_stencil end

function generate_kspace_stencil(
    recip_lattice::Mat3,
    kgrid_size::AbstractVector,
    kpoints::AbstractVector,
    ::FirstOrderKspaceStencil;
    atol=default_w90_kmesh_tol(),
)
    # find shells
    shells = search_shells(recip_lattice, kgrid_size, kpoints; atol)
    keep_shells = check_parallel(shells)
    shells = delete_shells(shells, keep_shells)

    keep_shells, bweights = compute_bweights(shells; atol)
    shells = delete_shells(shells, keep_shells)
    shells.bweights .= bweights

    check_completeness(shells; atol)
    # generate bvectors for each kpoint
    return sort_bvectors(shells; atol)
end

function generate_kspace_stencil(
    recip_lattice::Mat3, kgrid_size::AbstractVector, kpoints::AbstractVector; kwargs...
)
    return generate_kspace_stencil(
        recip_lattice, kgrid_size, kpoints, FirstOrderKspaceStencil(); kwargs...
    )
end

"""
    $(SIGNATURES)

Given bvector `b` connecting kpoints `ik` and `ikpb`, return the index of the
bvector `ib`.

This is a reverse search of bvector index if you only know the two kpoints
`ik` and `ikpb`, and the connecting displacement vector `G`, such that
`kpoints[ikpb] + G = kpoints[ik] + b`.

# Arguments
- `kpb_k`: length-`n_kpoints` vector, each element is a length-`n_bvectors`
    vector of k+b kpoint indices at kpoint `ik`
- `kpb_G`: length-`n_kpoints` vector, each element is a length-`n_bvectors`
    vector, then each element is a `Vec3` displacement vector for k+b bvectors
    at `ik`
- `ik`: integer, index of kpoint `k`
- `ikpb`: integer, index of kpoint `k+b`
- `G`: displacement vector from `k1` to `k2`, e.g. `Vec3{Int}`
"""
function index_bvector(
    kpb_k::AbstractVector,
    kpb_G::AbstractVector,
    ik::Integer,
    ikpb::Integer,
    G::AbstractVector,
)
    for (jb, (jk, jG)) in enumerate(zip(kpb_k[ik], kpb_G[ik]))
        if jk == ikpb && jG == G
            return jb
        end
    end
    return error("No neighbors found for ik = $(ik), ikpb = $(ikpb), G = $(G)")
end

"""
    $(SIGNATURES)

See also [`index_bvector`](@ref).
"""
function index_bvector(kstencil::KspaceStencil, ik, ikpb, G)
    return index_bvector(kstencil.kpb_k, kstencil.kpb_G, ik, ikpb, G)
end

function generate_kspace_stencil(
    recip_lattice::Mat3,
    kgrid_size::AbstractVector,
    kpoints::AbstractVector,
    ::CubicNearestKspaceStencil,
)
    dkx, dky, dkz = 1 ./ size(kgrid_size)
    T = eltype(recip_lattice)

    # only 6 nearest neighbors
    nbvecs = 6
    bvectors_frac = zeros(Vec3{T}, nbvecs)
    bvectors_frac[1] = Vec3{T}(dkx, 0, 0)
    bvectors_frac[2] = Vec3{T}(-dkx, 0, 0)
    bvectors_frac[3] = Vec3{T}(0, dky, 0)
    bvectors_frac[4] = Vec3{T}(0, -dky, 0)
    bvectors_frac[5] = Vec3{T}(0, 0, dkz)
    bvectors_frac[6] = Vec3{T}(0, 0, -dkz)

    # just a fake weight
    bweights = zeros(T, nbvecs)

    # generate bvectors for each kpoint, always the same across kpoints
    nkpts = length(kpoints)
    kpb_k = [zeros(Int, nbvecs) for _ in 1:nkpts]
    kpb_G = [zeros(Vec3{Int}, nbvecs) for _ in 1:nkpts]

    for (ik, kpt) in enumerate(kpoints)
        # use fractional coordinates for comparisons
        k_equiv, G_equiv = bvectors_to_kpb(bvectors_frac, kpt, kpoints)
        kpb_k[ik] = k_equiv
        kpb_G[ik] = G_equiv
    end

    bvectors = map(bvectors_frac) do b
        recip_lattice * b
    end
    return KspaceStencil(
        recip_lattice, kgrid_size, kpoints, bvectors, bweights, kpb_k, kpb_G
    )
end
