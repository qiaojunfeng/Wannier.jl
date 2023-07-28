using OffsetArrays: OffsetArray
using NearestNeighbors: KDTree, knn

export generate_Rspace, WignerSeitzRspace, MDRSRspace, BareRspace

"""
The Fourier-space frequencies for Wannier interpolation, also called
``\\mathbf{R}``-vectors.

# Fields

Since julia does not permit inheritance for fields, we need to ask the developers
to implement the following fields when introducing a new type of R-space domain,
see [`WignerSeitzRspace`](@ref) for an example.

- `lattice`: `3 * 3` matrix, each column is a lattice vector in Å unit
- `Rvectors`: length-`n_Rvectors` of `Vec3` for fractional (actually integers)
    coordinates w.r.t. lattice
"""
abstract type AbstractRspace end

n_Rvectors(Rspace::AbstractRspace) = length(Rspace.Rvectors)
real_lattice(Rspace::AbstractRspace) = Rspace.lattice

"""Index using `i` of Rvectors"""
function Base.getindex(Rspace::AbstractRspace, i::Integer)
    return Rspace.Rvectors[i]
end

Base.lastindex(Rspace::AbstractRspace) = lastindex(Rspace.Rvectors)
Base.length(Rspace::AbstractRspace) = length(Rspace.Rvectors)

function Base.iterate(Rspace::AbstractRspace, state=1)
    if state > length(Rspace.Rvectors)
        return nothing
    else
        return (Rspace.Rvectors[state], state + 1)
    end
end

function Base.show(io::IO, ::MIME"text/plain", Rspace::AbstractRspace)
    @printf(io, "R-space type  :  %s\n\n", nameof(typeof(Rspace)))
    show_lattice(io, Rspace.lattice)
    println(io)
    @printf(io, "n_Rvectors  =  %d", n_Rvectors(Rspace))
end

function Base.isapprox(a::AbstractRspace, b::AbstractRspace; kwargs...)
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
    $(TYPEDEF)

``\\mathbf{R}``-vectors generated using Wigner-Seitz cell.

# Fields
$(FIELDS)

!!! note

    The R vectors are sorted in the same order as `Wannier90`.
"""
struct WignerSeitzRspace{T<:Real} <: AbstractRspace
    """lattice, 3 * 3, each column is a lattice vector in Å unit"""
    lattice::Mat3{T}

    """R-vectors, length-`n_Rvectors` of `Vec3` for fractional (actually integers)
    coordinates w.r.t. lattice"""
    Rvectors::Vector{Vec3{Int}}

    """degeneracy of each Rvector, length-`n_Rvectors` vector.
    The weight of each Rvector = 1 / degeneracy."""
    n_Rdegens::Vector{Int}
end

function WignerSeitzRspace(
    lattice::AbstractMatrix, Rvectors::AbstractVector, n_Rdegens::AbstractVector
)
    T = eltype(lattice)
    return WignerSeitzRspace{T}(Mat3(lattice), Vector{Vec3{Int}}(Rvectors), n_Rdegens)
end

function WignerSeitzRspace(
    lattice::AbstractMatrix,
    Rgrid_size::Union{AbstractVector,Tuple};
    atol=default_w90_ws_distance_tol(),
    max_cell::Integer=default_w90_ws_search_size(),
)
    @assert length(Rgrid_size) == 3 "Rgrid_size should be a length-3 vector"
    # 1. generate a supercell where WFs live in
    supercell_wf, _ = make_supercell([Vec3(0, 0, 0)], [0:(r - 1) for r in Rgrid_size])
    # another supercell of the supercell_wf to find the Wigner-Seitz cell of the supercell_wf
    supercell, translations = make_supercell(
        supercell_wf, [((-max_cell):max_cell) * r for r in Rgrid_size]
    )
    # sort so z increases fastest, to make sure the Rvec order is the same as W90
    supercell = sort_points(supercell)
    # to cartesian coordinates
    supercell_cart = map(x -> lattice * x, supercell)
    # get translations of supercell_wf, only need unique points
    translations = unique(translations)
    translations_cart = map(x -> lattice * x, translations)

    # 2. KDTree to get the distance of supercell points to translations of lattice_wf
    kdtree = KDTree(translations_cart)
    # In priciple, need to calculate distances to all the supercell translations to
    # count degeneracies, this need a search of `size(translations_cart, 2)` neighbors.
    # usually we don't need such high degeneracies, so I only search for 8 neighbors.
    max_neighbors = min(8, length(translations_cart))
    idxs, dists = knn(kdtree, supercell_cart, max_neighbors, true)

    # 3. supercell_cart point which is closest to lattice_wf at origin is inside WS cell
    idx_origin = findfirst(isequal(Vec3(0, 0, 0)), translations)
    R_idxs = Vector{Int}()
    n_Rdegens = Vector{Int}()

    for (iR, R) in enumerate(supercell_cart)
        i = idxs[iR][1]
        d = dists[iR][1]
        if i != idx_origin
            # check again the distance, to reproduce W90's behavior
            if abs(d - norm(R)) >= atol
                continue
            end
        end
        push!(R_idxs, iR)
        degen = count(x -> abs(x - d) < atol, dists[iR])
        if degen > max_neighbors
            error("R-vector degeneracy is too large? $degen")
        end
        push!(n_Rdegens, degen)
    end
    # fractional coordinates
    Rvectors = supercell[R_idxs]

    return WignerSeitzRspace(lattice, Rvectors, n_Rdegens)
end

"""
    $(TYPEDEF)

The R-vectors and T-vectors for minimal-distance replica selection (MDRS) interpolation.

# Fields
$(FIELDS)

!!! note

    The R-vectors are sorted in the same order as wannier90.
"""
struct MDRSRspace{T<:Real} <: AbstractRspace
    """lattice, 3 * 3, each column is a lattice vector in Å unit"""
    lattice::Mat3{T}

    """R-vectors, length-`n_Rvectors` of `Vec3` for fractional (actually integers)
    coordinates w.r.t. lattice"""
    Rvectors::Vector{Vec3{Int}}

    """degeneracy of each Rvector, length-`n_Rvectors` vector.
    The weight of each Rvector = 1 / degeneracy."""
    n_Rdegens::Vector{Int}

    """translation T-vectors, fractional coordinates w.r.t lattice.
    Length-`n_Rvectors` vector, each element is a matrix of size
    `n_wannier * n_wannier`, then each element is a length-`n_Tdegens` vector
    of 3-vector for fractional coordinates"""
    Tvectors::Vector{Matrix{Vector{Vec3{Int}}}}

    """degeneracy of each T vector. Length-`n_Rvectors` vector, each element is
    a `n_wannier * n_wannier` matrix of integers"""
    n_Tdegens::Vector{Matrix{Int}}
end

function MDRSRspace(
    lattice::AbstractMatrix,
    Rvectors::AbstractVector,
    n_Rdegens::AbstractVector,
    Tvectors::AbstractVector,
    n_Tdegens::AbstractVector,
)
    T = eltype(lattice)
    return MDRSRspace{T}(Mat3(lattice), Rvectors, n_Rdegens, Tvectors, n_Tdegens)
end

function MDRSRspace(
    wsRspace::WignerSeitzRspace,
    Rgrid_size::Union{AbstractVector,Tuple},
    centers::AbstractVector;
    atol=default_w90_ws_distance_tol(),
    max_cell::Integer=default_w90_ws_search_size(),
)
    @assert length(Rgrid_size) == 3 "Rgrid_size should be a length-3 vector"
    nwann = length(centers)
    @assert nwann > 0 "centers is empty"

    nRvecs = n_Rvectors(wsRspace)
    lattice = wsRspace.lattice

    # 1. generate WS cell around origin to check WF |nR> is inside |m0> or not
    # increase max_cell by 1 in case WF center drifts away from the parallelepiped
    max_cell_1 = max_cell + 1
    # supercell of the supercell_wf to find the Wigner-Seitz cell of the supercell_wf,
    # supercell_wf is the cell where WFs live in
    supercell, translations = make_supercell(
        [Vec3(0, 0, 0)], [((-max_cell_1):max_cell_1) * r for r in Rgrid_size]
    )
    # to cartesian coordinates
    supercell_cart = map(x -> lattice * x, supercell)
    # get translations of supercell_wf, only need unique points
    translations = unique(translations)
    translations_cart = map(x -> lattice * x, translations)

    # 2. KDTree to get the distance of supercell points to translations of lattice_wf
    kdtree = KDTree(translations_cart)
    # usually we don't need such high degeneracies, so I only search for 8 neighbors
    max_neighbors = min(8, length(translations_cart))
    idx_origin = findfirst(isequal(Vec3(0, 0, 0)), translations)
    # save all translations and degeneracies
    Tvectors = [Matrix{Vector{Vec3{Int}}}(undef, nwann, nwann) for _ in 1:nRvecs]
    Tdegens = [zeros(Int, nwann, nwann) for _ in 1:nRvecs]

    @inbounds @views for iR in 1:nRvecs
        Tvectors_iR = Tvectors[iR]
        Tdegens_iR = Tdegens[iR]
        R = wsRspace.Rvectors[iR]
        for m in 1:nwann
            for n in 1:nwann
                # translation vector of |nR> WF center relative to |m0> WF center
                Tᶠ = centers[n] + R .- centers[m]
                # to cartesian
                Tᶜ = map(x -> x .+ lattice * Tᶠ, supercell_cart)
                # get distances
                idxs, dists = knn(kdtree, Tᶜ, max_neighbors, true)
                # collect T vectors
                T_idxs = Vector{Int}()
                for iT in 1:length(Tᶜ)
                    i = idxs[iT][1]
                    d = dists[iT][1]
                    if i != idx_origin
                        # check again the distance, to reproduce W90's behavior
                        if idx_origin ∉ idxs[iT]
                            continue
                        end
                        j = findfirst(idxs[iT] .== idx_origin)
                        d0 = dists[iT][j]
                        if abs(d - d0) >= atol
                            continue
                        end
                    end
                    push!(T_idxs, iT)
                end
                degen = length(T_idxs)
                if degen > max_neighbors
                    error("degeneracy of T-vectors is too large? $degen")
                end
                # fractional coordinates
                Tvectors_iR[m, n] = map(i -> supercell[i], T_idxs)
                Tdegens_iR[m, n] = degen
            end
        end
    end
    return MDRSRspace(
        wsRspace.lattice, wsRspace.Rvectors, wsRspace.n_Rdegens, Tvectors, Tdegens
    )
end

function MDRSRspace(
    lattice::AbstractMatrix,
    Rgrid_size::Union{AbstractVector,Tuple},
    centers::AbstractVector;
    kwargs...,
)
    wsRspace = WignerSeitzRspace(lattice, Rgrid_size; kwargs...)
    return MDRSRspace(wsRspace, Rgrid_size, centers; kwargs...)
end

"""
    $(SIGNATURES)

Generate R-space domain for Wannier interpolation.

# Arguments
- `lattice`: columns are lattice vectors
- `Rgrid_size`: number of FFT grid points in each direction, actually determined
    by (and is equal to) `kgrid_size`
- `centers`: length-`n_wannier` vector, each element is a WF center in
    fractional coordinates

# Keyword arguments
- `atol`: toerance for checking degeneracy
- `max_cell`: number of neighboring cells to be searched

# Return
- `Rdomain`: a [`WignerSeitzRspace`](@ref) or [`MDRSRspace`](@ref)

!!! note

    To reproduce wannier90's behavior
    - `atol` should be equal to wannier90's input parameter `ws_distance_tol`
    - `max_cell` should be equal to wannier90's input parameter `ws_search_size`
"""
function generate_Rspace end

function generate_Rspace(T::Type{<:AbstractRspace}, args...; kwargs...)
    return T(args...; kwargs...)
end

function generate_Rspace(model::Model; MDRS::Bool=true, kwargs...)
    if MDRS
        # from Cartesian to fractional
        inv_lattice = inv(model.lattice)
        r = map(center(model)) do c
            inv_lattice * c
        end
        Rspace = MDRSRspace(model.lattice, model.kgrid_size, r; kwargs...)
    else
        Rspace = WignerSeitzRspace(model.lattice, model.kgrid_size; kwargs...)
    end
    return Rspace
end

"""Type for the `xyz_iR` mapping using OffsetArray"""
const RvectorIndexMapping = OffsetArray{Int,3}

"""
    $(TYPEDEF)

A minimalistic R-space domain, with only R-vectors themselves.

Contrary to [`WignerSeitzRspace`](@ref) and [`MDRSRspace`](@ref), this domain
does not contain any info on R-vector degeneracies or translation T-vectors.
Therefore, it is not possible to do the forward Fourier transform from k-space
operators to R-space representation (in comparison, the previous two domains
contain such info for the forward Fourier transform). However, we can absorb
the degeneracies and T-vectors into the R-space operators themselves during
forward Fourier transform, and the backward Fourier transforms are just simple
sum over ``\\exp(i \\mathbf{k} \\cdot \\mathbf{R})``, reducing the computational
cost. This is why we name it "bare" R-space domain. Using this domain, the
R-space operators are straightforward tight-binding models, expressed as
Fourier series.

# Fields
$(FIELDS)
"""
struct BareRspace{T<:Real} <: AbstractRspace
    """lattice, 3 * 3, each column is a lattice vector in Å unit"""
    lattice::Mat3{T}

    """R-vectors, length-`n_Rvectors` of `Vec3` for fractional (actually integers)
    coordinates w.r.t. lattice"""
    Rvectors::Vector{Vec3{Int}}

    """the mapping from `[i, j, k]` index to `iR` index, to allow indexing
    the struct using the x,y,z component of R-vectors. This can be auto
    constructed using [`build_mapping_xyz_iR`](@ref)."""
    xyz_iR::RvectorIndexMapping
end

"""
    $(SIGNATURES)

Build the mapping from R-vector itself to its index `iR`, such that
`mapping[Rvectors[iR]...] == iR`. Missing indices are filled with `0`.
"""
function build_mapping_xyz_iR(Rvectors::AbstractVector)
    Rx = [R[1] for R in Rvectors]
    Ry = [R[2] for R in Rvectors]
    Rz = [R[3] for R in Rvectors]
    minRx, maxRx = extrema(Rx)
    minRy, maxRy = extrema(Ry)
    minRz, maxRz = extrema(Rz)
    nRx = maxRx - minRx + 1
    nRy = maxRy - minRy + 1
    nRz = maxRz - minRz + 1
    # fill the mapping with 0 to indicate that index is not valid,
    # cannot use NaN or Inf because they are not floats
    mapping = zeros(Int, nRx, nRy, nRz)
    # to offset array
    mapping = RvectorIndexMapping(mapping, minRx:maxRx, minRy:maxRy, minRz:maxRz)
    for (iR, Rvec) in enumerate(Rvectors)
        mapping[Rvec...] = iR
    end
    # to StaticArray so that the indices can not be modified
    mapping = SArray{Tuple{size(mapping)...}}(mapping)
    mapping = RvectorIndexMapping(mapping, minRx:maxRx, minRy:maxRy, minRz:maxRz)
    return mapping
end

function BareRspace{T}(lattice::AbstractMatrix, Rvectors::AbstractVector) where {T}
    xyz_iR = build_mapping_xyz_iR(Rvectors)
    return BareRspace{T}(Mat3(lattice), Vector{Vec3{Int}}(Rvectors), xyz_iR)
end

function BareRspace(lattice::AbstractMatrix, Rvectors::AbstractVector)
    T = eltype(lattice)
    return BareRspace{T}(lattice, Rvectors)
end
