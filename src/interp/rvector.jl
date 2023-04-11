using LinearAlgebra
using NearestNeighbors: KDTree, knn

export get_Rvectors_ws, get_Rvectors_mdrs

"""
    struct RVectors

The R vectors for interpolation.

# Fields
- `lattice`: columns are lattice vectors
- `grid`: number of FFT grid points in each direction, actually equal to `kgrid`
- `R`: `3 * n_rvecs`, R vectors in fractional coordinates w.r.t. lattice
- `N`: `n_rvecs`, degeneracy of each R vector

!!! note

    The R vectors are sorted in the same order as `Wannier90`.
"""
struct RVectors{T<:Real}
    # lattice, 3 * 3, Å unit
    # each column is a lattice vector
    lattice::Mat3{T}

    # grid size along x,y,z, actually equal to kgrid
    grid::Vec3{Int}

    # R vectors, 3 * n_rvecs
    # fractional (actually integers) coordinates w.r.t lattice
    R::Matrix{Int}

    # degeneracy of each Rvec, n_rvecs
    # weight = 1 / degeneracy
    N::Vector{Int}
end

function Base.getproperty(x::RVectors, sym::Symbol)
    if sym == :n_rvecs
        return size(x.R, 2)
    else
        # fallback to getfield
        getfield(x, sym)
    end
end

"""
    RVectors(lattice, grid, R, N)

Constructor for `RVectors`.

Auto transform `lattice` and `grid` to `Mat3` and `Vec3`, respectively.

# Arguments
- `lattice`: columns are lattice vectors
- `grid`: number of FFT grid points in each direction, actually equal to `kgrid`
- `R`: `3 * n_rvecs`, R vectors in fractional coordinates w.r.t. lattice
- `N`: `n_rvecs`, degeneracy of each R vector
"""
function RVectors(
    lattice::AbstractMatrix,
    grid::AbstractVector{T},
    R::AbstractMatrix{T},
    N::AbstractVector{T},
) where {T<:Integer}
    return RVectors(Mat3(lattice), Vec3(grid), R, N)
end

function Base.show(io::IO, Rvectors::RVectors)
    @printf(io, "lattice: Å\n")
    for i in 1:3
        @printf(io, "  a%d: %8.5f %8.5f %8.5f\n", i, Rvectors.lattice[:, i]...)
    end
    @printf(io, "grid    = %d %d %d\n", Rvectors.grid...)
    @printf(io, "n_rvecs = %d", Rvectors.n_rvecs)
end

"""
    check_weights(R::RVectors)

Sanity check for the degeneracies of R vectors.
"""
function check_weights(R::RVectors)
    if sum(1 ./ R.N) ≉ prod(R.grid)
        error("weights do not sum to 1")
    end
end

"""
    get_Rvectors_ws(lattice, rgrid; atol=1e-5, max_cell=3)

Generate R vectors for Wigner-Seitz interpolation.

# Arguments
- `lattice`: columns are lattice vectors
- `rgrid`: number of FFT grid points in each direction, actually equal to `kgrid`

# Keyword arguments
- `atol`: toerance for checking degeneracy,
    equivalent to `Wannier90` input parameter `ws_distance_tol`
- `max_cell`: number of neighboring cells to be searched,
    equivalent to `Wannier90` input parameter `ws_search_size`
"""
function get_Rvectors_ws(
    lattice::AbstractMatrix{T}, rgrid::AbstractVector{R}; atol::T=1e-5, max_cell::Int=3
) where {T<:Real,R<:Integer}
    # 1. Generate a supercell where WFs live in
    supercell_wf, _ = make_supercell(zeros(Int, 3, 1), [0:(r - 1) for r in rgrid])
    # another supercell of the supercell_wf to find the Wigner Seitz cell of the supercell_wf
    supercell, translations = make_supercell(
        supercell_wf, [((-max_cell):max_cell) * r for r in rgrid]
    )
    # sort so z increases fastest, to make sure the Rvec order is the same as W90
    supercell = sort_kpoints(supercell)
    # to cartesian coordinates
    supercell_cart = lattice * supercell
    # get translations of supercell_wf, only need unique points
    translations = unique(translations; dims=2)
    translations_cart = lattice * translations

    # 2. KDTree to get the distance of supercell points to translations of lattice_wf
    kdtree = KDTree(translations_cart)
    # in priciple, need to calculate distances to all the supercell translations to
    # count degeneracies, this need a search of `size(translations_cart, 2)` neighbors.
    # usually we don't need such high degeneracies, so I only search for 8 neighbors.
    max_neighbors = min(8, size(translations_cart, 2))
    idxs, dists = knn(kdtree, supercell_cart, max_neighbors, true)

    # 3. supercell_cart point which is closest to lattice_wf at origin is inside WS cell
    idx_origin = findvector(==, [0, 0, 0], translations)
    R_idxs = Vector{Int}()
    R_degen = Vector{Int}()
    for ir in axes(supercell_cart, 2)
        i = idxs[ir][1]
        d = dists[ir][1]
        if i != idx_origin
            # check again the distance, to reproduce W90's behavior
            d0 = norm(supercell_cart[:, ir])
            if abs(d - d0) >= atol
                continue
            end
        end
        push!(R_idxs, ir)
        degen = count(x -> abs(x - d) < atol, dists[ir])
        if degen > max_neighbors
            error("degeneracy is too large? $degen")
        end
        push!(R_degen, degen)
    end
    # fractional coordinates
    R_vecs = supercell[:, R_idxs]

    Rvecs = RVectors(lattice, rgrid, R_vecs, R_degen)
    check_weights(Rvecs)

    return Rvecs
end

"""
    struct RVectorsMDRS

The R vectors for MDRS interpolation.

# Fields
- `Rvectors`: `Rvectors` for Wigner-Seitz interpolation
## For MDRSv1
- `T`: `n_wann * n_wann * n_rvecs`, translation vectors w.r.t to lattice for MDRSv1
- `Nᵀ`: `n_wann * n_wann * n_rvecs`, degeneracy of each `T` vector for MDRSv1
## For MDRSv2
- `R̃vectors`: `RVectors` containing expanded set of `R + T` vectors used in MDRSv2
- `R̃_RT`: mapping of `R̃vectors` to `R` and `T` vectors

!!! note

    The R vectors are sorted in the same order as `Wannier90`.
"""
struct RVectorsMDRS{U<:Real}
    # R vectors of the Wigner-Seitz interplation
    Rvectors::RVectors{U}

    # For MDRS v1

    # translation vectors, fractional coordinates w.r.t lattice
    # internal matrix: 3 * n_degen, external array: n_wann * n_wann * n_rvecs
    T::Array{Matrix{Int},3}

    # degeneracy of each T vector, n_wann * n_wann * n_rvecs
    Nᵀ::Array{Int,3}

    # For MDRS v2

    # R̃ vectors (expanded set for R+T), 3 * n_r̃vecs
    R̃vectors::RVectors{U}

    # mapping of R̃vectors to R and T vectors,
    # length = n_r̃vecs, each element is a Vector whose element are [ir, m, n, it],
    # where ir is the index of R vector, (m, n, it) specifies the T vector for
    # WFs |0m> and |Rn>
    R̃_RT::Vector{Vector{SVector{4,Int}}}
end

"""
    RVectorsMDRS(Rvectors, T, Nᵀ)

A friendly constructor for `RVectorsMDRS`.

The remaining fields `R̃vectors` and `R̃_RT` are only used in MDRSv2,
the are calculated automatically based on the input arguments.

# Arguments
- `Rvectors`: `RVectors` for Wigner-Seitz interpolation
- `T`: `n_wann * n_wann * n_rvecs`, translation vectors w.r.t to lattice for MDRSv1
- `Nᵀ`: `n_wann * n_wann * n_rvecs`, degeneracy of each `T` vector for MDRSv1
"""
function RVectorsMDRS(Rvectors::RVectors, T::Array{Matrix{Int},3}, Nᵀ::Array{Int,3})
    # expanded R vectors
    R̃ = Vector{Vec3{Int}}()
    # mapping
    R̃_RT = Vector{Vector{SVector{4,Int}}}()
    # generate expanded R̃ vectors, which contains all the R+T
    for ir in 1:(Rvectors.n_rvecs)
        for n in axes(T, 2)
            for m in axes(T, 1)
                for it in 1:Nᵀ[m, n, ir]
                    RT = Rvectors.R[:, ir] + T[m, n, ir][:, it]
                    i = findfirst(x -> x == RT, R̃)
                    if isnothing(i)
                        push!(R̃, RT)
                        push!(R̃_RT, [[ir, m, n, it]])
                    else
                        push!(R̃_RT[i], [ir, m, n, it])
                    end
                end
            end
        end
    end
    # I will multiply R degeneracy when doing fourier k -> R̃, contrary to WS interpolation
    # where R degeneracy is multiplied during inv fourier R -> k. Si here Ñ is 1.
    Ñ = zeros(Int, length(R̃))
    # to matrix
    R̃mat = zeros(Int, 3, length(R̃))
    for ir̃ in axes(R̃mat, 2)
        R̃mat[:, ir̃] = R̃[ir̃]
    end
    R̃vectors = RVectors(Rvectors.lattice, Rvectors.grid, R̃mat, Ñ)
    return RVectorsMDRS(Rvectors, T, Nᵀ, R̃vectors, R̃_RT)
end

function Base.getproperty(x::RVectorsMDRS, sym::Symbol)
    if sym ∈ fieldnames(RVectors)
        return getfield(x.Rvectors, sym)
    elseif sym == :n_rvecs
        return getproperty(x.Rvectors, sym)
    elseif sym == :n_r̃vecs
        return x.R̃vectors.n_rvecs
    else
        # fallback to getfield
        getfield(x, sym)
    end
end

function Base.show(io::IO, Rvectors::RVectorsMDRS)
    @printf(io, "lattice: Å\n")
    for i in 1:3
        @printf(io, "  a%d: %8.5f %8.5f %8.5f\n", i, Rvectors.lattice[:, i]...)
    end
    @printf(io, "grid    = %d %d %d\n", Rvectors.grid...)
    @printf(io, "n_rvecs = %d\n", Rvectors.n_rvecs)
    return print(io, "using MDRS interpolation")
end

"""
    get_Rvectors_mdrs(lattice, rgrid, centers; atol=1e-5, max_cell=3)

Generate R vectors for MDRS interpolation (both v1 and v2).

# Arguments
- `lattice`: columns are lattice vectors
- `rgrid`: number of FFT grid points in each direction, actually equal to `kgrid`
- `centers`: `3 * n_wann`, WF centers in fractional coordinates

# Keyword arguments
- `atol`: toerance for checking degeneracy,
    equivalent to `Wannier90` input parameter `ws_distance_tol`
- `max_cell`: number of neighboring cells to be searched,
    equivalent to `Wannier90` input parameter `ws_search_size`
"""
function get_Rvectors_mdrs(
    lattice::AbstractMatrix{T},
    rgrid::AbstractVector{Int},
    centers::AbstractMatrix{T};
    atol::T=1e-5,
    max_cell::Int=3,
) where {T<:Real}
    n_wann = size(centers, 2)
    Rvec = get_Rvectors_ws(lattice, rgrid; atol=atol, max_cell=max_cell)
    n_rvecs = Rvec.n_rvecs

    # 1. generate WS cell around origin to check WF |nR> is inside |m0> or not
    # increase max_cell by 1 in case WF center drifts away from the parallelepiped
    max_cell1 = max_cell + 1
    # supercell of the supercell_wf to find the Wigner Seitz cell of the supercell_wf,
    # supercell_wf is the cell where WFs live in
    supercell, translations = make_supercell(
        zeros(Int, 3, 1), [((-max_cell1):max_cell1) * r for r in rgrid]
    )
    # to cartesian coordinates
    supercell_cart = lattice * supercell
    # get translations of supercell_wf, only need unique points
    translations = unique(translations; dims=2)
    translations_cart = lattice * translations

    # 2. KDTree to get the distance of supercell points to translations of lattice_wf
    kdtree = KDTree(translations_cart)
    # usually we don't need such high degeneracies, so I only search for 8 neighbors
    max_neighbors = min(8, size(translations_cart, 2))
    idx_origin = findvector(==, [0, 0, 0], translations)

    # save all translations and degeneracies
    T_vecs = Array{Matrix{Int}}(undef, n_wann, n_wann, n_rvecs)
    T_degen = zeros(Int, n_wann, n_wann, n_rvecs)

    for ir in 1:n_rvecs
        for m in 1:n_wann
            for n in 1:n_wann
                # translation vector of |nR> WFC relative to |m0> WFC
                Tᶠ = centers[:, n] + Rvec.R[:, ir] - centers[:, m]
                # to cartesian
                Tᶜ = supercell_cart .+ lattice * Tᶠ
                # get distances
                idxs, dists = knn(kdtree, Tᶜ, max_neighbors, true)
                # collect T vectors
                T_idxs = Vector{Int}()
                for it in axes(Tᶜ, 2)
                    i = idxs[it][1]
                    d = dists[it][1]
                    if i != idx_origin
                        # check again the distance, to reproduce W90's behavior
                        if idx_origin ∉ idxs[it]
                            continue
                        end
                        j = findfirst(idxs[it] .== idx_origin)
                        d0 = dists[it][j]
                        if abs(d - d0) >= atol
                            continue
                        end
                    end
                    push!(T_idxs, it)
                end
                degen = length(T_idxs)
                if degen > max_neighbors
                    error("degeneracy of T vectors is too large? $degen")
                end
                # fractional coordinates
                T_vecs[m, n, ir] = supercell[:, T_idxs]
                T_degen[m, n, ir] = degen
            end
        end
    end

    return RVectorsMDRS(Rvec, T_vecs, T_degen)
end
