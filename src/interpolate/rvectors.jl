using LinearAlgebra
import NearestNeighbors as NN

@doc raw"""
The R vectors for interpolations, sorted in the same order as the W90 nnkp file.
"""
struct RVectors{T<:Real}
    # lattice, 3 * 3, Å unit
    # each column is a lattice vector
    lattice::Mat3{T}

    # grid size along x,y,z, actually equal to kgrid
    grid::Vec3{Int}

    # R vectors, Cartesian! coordinates, Å unit, 3 * n_rvecs
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

function RVectors(
    lattice::AbstractMatrix,
    grid::AbstractVector{T},
    R::AbstractMatrix{T},
    N::AbstractVector{T},
) where {T<:Integer}
    return RVectors(Mat3(lattice), Vec3(grid), R, N)
end

function check_weights(R::RVectors)
    if sum(1 ./ R.N) ≉ prod(R.grid)
        error("weights do not sum to 1")
    end
end

"""
atol: equivalent to `ws_distance_tol` in wannier90.
"""
function get_Rvectors_ws(
    lattice::AbstractMatrix{T}, rgrid::AbstractVector{Int}; atol::T=1e-5, max_cell::Int=3
) where {T<:Real}
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
    kdtree = NN.KDTree(translations_cart)
    # need to calculate distances to all the supercell translations to count degeneracies
    idxs, dists = NN.knn(kdtree, supercell_cart, size(translations_cart, 2), true)

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
        push!(R_degen, degen)
    end
    # fractional coordinates
    R_vecs = supercell[:, R_idxs]

    R = RVectors(lattice, rgrid, R_vecs, R_degen)
    check_weights(R)

    return R
end

function get_Rvectors_mdrs(lattice::AbstractMatrix, kgrid::AbstractVector{Int}) end
