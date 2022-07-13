using Printf: @printf
import LinearAlgebra as LA
import NearestNeighbors as NN

@doc raw"""
Store shells of bvectors.
"""
struct BVectorShells{T<:Real}
    # reciprocal lattice, 3 * 3
    recip_lattice::Mat3{T}

    # kpoints array, fractional coordinates, n_kpts of Vec3
    kpoints::Matrix{T}

    # bvectors of each shell, Cartesian coordinates, n_shells of 3 * multiplicity
    bvectors::Vector{Matrix{T}}

    # weight of each shell, n_shells
    weights::Vector{T}

    # multiplicity of each shell, num_shells
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

@doc raw"""
The bvectors for each kpoint, sorted the same as the Wannier90 nnkp order.
"""
struct BVectors{T<:Real}
    # reciprocal lattice, 3 * 3
    recip_lattice::Mat3{T}

    # kpoints array, fractional coordinates, n_kpts of Vec3
    kpoints::Matrix{T}

    # bvectors, Cartesian coordinates, 3 * n_bvecs
    bvectors::Matrix{T}

    # weight of each bvec, n_bvecs
    weights::Vector{T}

    # k+b vectors, k -> k + b (index of equivalent kpt in the 1st BZ), n_bvecs * n_kpts
    kpb_k::Matrix{Int}

    # displacements between k + b and k + b wrapped around into the recip_cell,
    # 3 * n_bvecs * n_kpts, where 3 is [b_x, b_y, b_z]
    kpb_b::Array{Int,3}
end

function Base.getproperty(x::BVectors, sym::Symbol)
    if sym == :n_kpts
        return length(x.kpoints)
    elseif sym == :n_bvecs
        return size(x.bvectors, 2)
    else
        # fallback to getfield
        getfield(x, sym)
    end
end

@doc raw"""
Make a supercell of kpoints by translating it along 3 directions.
Input and returned kpoints are in fractional coordinates.
"""
function make_supercell(kpoints::Matrix{T}, supercell_size::Int=5) where {T<:Real}
    size(kpoints, 1) ≉ 3 && error("kpoints must be 3 * n_kpts")
    n_kpts = size(kpoints, 2)

    n_cells = (2 * supercell_size + 1)^3

    supercell = Matrix{T}(undef, 3, n_cells * n_kpts)
    translations = Matrix{Int}(undef, 3, n_cells * n_kpts)

    counter = 1
    # Note: the order of which index increases the fastest is crucial,
    # it will determine the order of bvecotors.
    for ix in (-supercell_size):supercell_size
        for iy in (-supercell_size):supercell_size
            for iz in (-supercell_size):supercell_size
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
Check if the columns of matrix A and columns of matrix B are parallel.
"""
function are_parallel(A::Matrix{T}, B::Matrix{T}; atol::T=1e-5) where {T<:Real}
    n_dim = size(A, 1)
    if n_dim != 3 || size(B, 1) != n_dim
        error("only support 3-vectors")
    end

    ncol_A = size(A, 2)
    ncol_B = size(B, 2)
    result = fill(false, ncol_A, ncol_B)

    for (i, col1) in enumerate(eachcol(A))
        for (j, col2) in enumerate(eachcol(B))
            p = LA.cross(col1, col2)
            if all(isapprox.(0, p; atol=atol))
                result[i, j] = true
            end
        end
    end

    return result
end

@doc raw"""
Try to guess bvector weights from MV1997 Eq. (B1).

input bvectors are overcomplete vectors found during shell search.
return: bvectors and weights satisfing B1 condition.
"""
function get_weights(bvectors::Vector{Matrix{T}}; b1_atol::T=1e-6) where {T<:Real}
    n_shells = length(bvectors)
    n_shells == 0 && error("empty bvectors?")
    size(bvectors[1], 1) != 3 && error("only support 3-vectors")

    # only compare the upper triangular part of bvec * bvec', 6 elements
    B = zeros(T, 6, n_shells)

    # return the upper triangular part of a matrix as a vector
    triu2vec(m) = m[LA.triu!(trues(size(m)), 0)]
    # triu2vec(I) = [1 0 1 0 0 1]
    triu_I = triu2vec(LA.diagm([1, 1, 1]))

    W = zeros(n_shells)

    # sigular value tolerance
    σ_atol = 1e-5

    keep_shells = zeros(Int, 0)
    ishell = 1
    while ishell <= n_shells
        push!(keep_shells, ishell)
        B[:, ishell] = triu2vec(bvectors[ishell] * bvectors[ishell]')
        # Solve equation B * W = triu_I
        # size(B) = (6, ishell), W is diagonal matrix of size ishell
        # B = U * S * V' -> W = V * S^-1 * U' * triu_I
        U, S, V = LA.svd(B[:, keep_shells])
        @debug "S" ishell S = S' keep_shells = keep_shells'
        if all(S .> σ_atol)
            W[keep_shells] = V * LA.Diagonal(S)^-1 * U' * triu_I
            BW = B[:, keep_shells] * W[keep_shells]
            @debug "BW" ishell BW = BW'
            if isapprox(BW, triu_I; atol=b1_atol)
                break
            end
        else
            pop!(keep_shells)
        end
        ishell += 1
    end
    if ishell == n_shells + 1
        error("not enough shells to satisfy B1 condition")
    end
    n_shells = length(keep_shells)

    # resize
    weights = W[keep_shells]
    new_bvectors = bvectors[keep_shells]

    return new_bvectors, weights
end

@doc raw"""
Search bvector shells satisfing B1 condition.

kpoints: fractional coordinates
recip_lattice: each column is a reciprocal lattice vector
atol: precision to select a shell (equal distance to the point)
"""
function search_shells(
    kpoints::Matrix{T}, recip_lattice::Mat3{T}; atol::T=1e-5
) where {T<:Real}
    # Usually these "magic" numbers work well for normal recip_lattice.
    # Number of nearest-neighbors to be returned
    max_neighbors = 100
    # Max number of stencils in one shell
    max_multiplicity = 40
    # Max number of nearest-neighbor shells
    max_shells = 5
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
    ineigh = 2
    ishell = 1
    while (ineigh <= max_neighbors) && (ishell <= max_shells)
        # use the 1st kpoint to find bvector shells & weights
        eqdist_idxs = findall(x -> isapprox(x, dists[ineigh]; atol=atol), dists)
        multi = length(eqdist_idxs)
        if multi >= max_multiplicity
            # skip large multiplicity shells
            ineigh += multi
            break
        end
        push!(shells, idxs[eqdist_idxs])
        ineigh += multi
        ishell += 1
    end

    # 4. Get Cartesian coordinates vectors
    n_shells = length(shells)
    bvectors = Vector{Matrix{T}}(undef, n_shells)
    for ishell in 1:n_shells
        kpb_cart = supercell_cart[:, shells[ishell]]
        bvectors[ishell] = kpb_cart .- kpt_orig
    end
    @debug "Found bvector shells" bvectors

    # 5. Remove shells with parallel bvectors
    keep_shells = collect(1:n_shells)
    for ishell in 2:n_shells
        has_parallel = false
        for jshell in 1:(ishell - 1)
            if !(jshell in keep_shells)
                continue
            end
            if has_parallel
                break
            end
            p = are_parallel(bvectors[jshell], bvectors[ishell])
            # @debug "are_parallel($jshell, $ishell)" p
            if any(p)
                has_parallel = true
                break
            end
        end
        if has_parallel
            splice!(keep_shells, ishell)
        end
    end
    n_shells = length(keep_shells)
    bvectors = bvectors[keep_shells]

    # @debug "After check_parallel bvector shells" bvectors

    # 6. Calculate weights to satisfy B1 condition
    bvectors, weights = get_weights(bvectors)
    n_shells = length(bvectors)

    for i in 1:n_shells
        @printf("BVector shell %3d    weight = %8.5f\n", i, weights[i])
        vecs = bvectors[i]
        for ib in 1:size(vecs, 2)
            @printf("  %3d    %10.5f %10.5f %10.5f\n", ib, vecs[:, ib]...)
        end
    end

    return BVectorShells(recip_lattice, kpoints, bvectors, weights)
end

@doc raw"""
Check completeness (B1 condition) of bvectors.
"""
function check_b1(shells::BVectorShells{T}, b1_atol::T=1e-6) where {T<:Real}
    M = zeros(T, 3, 3)

    for ish in 1:(shells.n_shells)
        bvec = shells.bvectors[ish]
        M += shells.weights[ish] * bvec * bvec'
    end

    @debug "Bvector sum" M
    if !isapprox(M, LA.I; atol=b1_atol)
        error("B1 condition is not satisfied")
    end

    println("Finite difference condition satisfied")
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
    for ishell in 1:(shells.n_shells)
        multi = shells.multiplicities[ishell]
        bvecs[:, counter:(counter + multi - 1)] = shells.bvectors[ishell]
        bvecs_weight[counter:(counter + multi - 1)] .= shells.weights[ishell]
        counter += multi
    end

    return bvecs, bvecs_weight
end

@doc raw"""
Sort bvectors in shells for each kpoints, to be consistent with wannier90

Wannier90 use different order of bvectors for each kpoint,
in principle, this is not needed. However, the Mmn file is
written in such order, so I need to sort bvectors and calculate
weights, since nnkp file has no section of weights.
"""
function sort_bvectors(shells::BVectorShells{T}) where {T<:Real}
    # Firstly, sort by length of bvectors: nearest k+b goes first
    # Secondly, sort by supercell: k+b which is in the same cell as k goes first
    # Thirdly, sort by the index (in the supercell array) of supercell:
    #     the larger-index one goes first.
    #     Note this also relies on the order of supercell_cart, see make_supercell()
    # Fourthly, sort by k+b index: the smaller-index k+b point goes first

    kpoints = shells.kpoints
    recip_lattice = shells.recip_lattice

    n_kpts = size(kpoints, 2)

    inv_recip = inv(recip_lattice)
    kpts_cart = recip_lattice * kpoints

    # To sort bvectors for each kpoints, I need to
    # calculate distances of supercells to original cell
    # Just pass a fake kpoint, I only need the cell translations.
    _, supercell_idx = make_supercell(zeros(T, 3, 1))

    """Equivalent to periodic image?"""
    function isequiv(v1, v2; atol=1e-5)
        d = v1 - v2
        d -= round.(d)

        return all(isapprox.(d, 0; atol=atol))
    end

    """Find vector in the columns of a matrix"""
    function findvector(predicate::Function, v, M)
        for (i, col) in enumerate(eachcol(M))
            predicate(v, col) && return i
        end
        return error("$v not found in array!")
    end

    """
    kb1, kb2, k: cartesian coordinates
    """
    function bvec_isless(kb1, kb2, k; atol=1e-5)
        # Float64 comparsion is tricky, esp. for norm

        if isapprox(kb1, kb2; atol=atol)
            return true
        end

        normkb1 = LA.norm(kb1 - k)
        normkb2 = LA.norm(kb2 - k)
        if abs(normkb1 - normkb2) < atol
            # using outter scope variable: inv_recip, kpoints, kpts_cart
            # bring back to fractional coord, easier to compare
            kb1_frac = inv_recip * kb1
            kb2_frac = inv_recip * kb2
            kb1_equiv_idx = findvector(isequiv, kb1_frac, kpoints)
            kb2_equiv_idx = findvector(isequiv, kb2_frac, kpoints)
            celldisp_kb1 = kb1 - kpts_cart[:, kb1_equiv_idx]
            celldisp_kb2 = kb2 - kpts_cart[:, kb2_equiv_idx]
            # @debug "bvec_isless" celldisp_kb1, celldisp_kb2
            normkb1 = LA.norm(celldisp_kb1)
            normkb2 = LA.norm(celldisp_kb2)
            # @debug "bvec_isless" normkb1, normkb2
            if abs(normkb1 - normkb2) < atol
                # using outter scope variable: supercell_idx
                celldisp_kb1_frac = round.(Int, inv_recip * celldisp_kb1)
                celldisp_kb2_frac = round.(Int, inv_recip * celldisp_kb2)
                idx1 = findvector(==, celldisp_kb1_frac, supercell_idx)
                idx2 = findvector(==, celldisp_kb2_frac, supercell_idx)
                # @debug "bvec_isless", idx1, idx2
                if idx1 > idx2
                    return true
                elseif idx1 < idx2
                    return false
                else
                    if kb1_equiv_idx < kb2_equiv_idx
                        return true
                    elseif kb1_equiv_idx > kb2_equiv_idx
                        return false
                    else
                        error("Comparing equivalent k+b points, this should not happen!")
                    end
                end
            elseif normkb1 < normkb2
                return true
            else # normkb1 > normkb2
                return false
            end
        elseif normkb1 < normkb2
            return true
        else # normkb1 > normkb2
            return false
        end
    end

    # Test: check the order of b vectors is the same as wannier90
    # k = kpts_cart[:, 1]
    # kb1 = recip_lattice * (kpoints[:, 9] + [0 0 0]')
    # kb2 = recip_lattice * (kpoints[:, 65] + [0 0 0]')
    # @debug "islesss" bvec_isless(kb1, kb2, k)
    # throw(ErrorException)

    bvecs, bvecs_weight = flatten_shells(shells)
    n_bvecs = size(bvecs, 2)

    # find k+b indexes
    kpb_k = zeros(Int, n_bvecs, n_kpts)
    kpb_weight = zeros(T, n_bvecs, n_kpts)
    kpb_b = zeros(Int, 3, n_bvecs, n_kpts)
    sorted_idx = zeros(Int, n_bvecs)
    # I have to construct a Vector of Vector such that sortperm could work
    ikpb = Vector{Vec3{T}}(undef, n_bvecs)

    for ik in 1:n_kpts
        k_cart = kpts_cart[:, ik]

        for ib in 1:n_bvecs
            ikpb[ib] = k_cart .+ bvecs[:, ib]
        end

        sorted_idx[:] = sortperm(ikpb; lt=(x, y) -> bvec_isless(x, y, k_cart))

        for ib in 1:n_bvecs
            ib_sorted = sorted_idx[ib]
            kpb_weight[ib, ik] = bvecs_weight[ib_sorted]

            ikpb_frac = inv_recip * ikpb[ib_sorted]
            ikpb_equiv_idx = findvector(isequiv, ikpb_frac, kpoints)
            kpb_k[ib, ik] = ikpb_equiv_idx

            kpb_b[:, ib, ik] = round.(Int, ikpb_frac - kpoints[:, ikpb_equiv_idx])
        end
    end

    # kpb_weight is redundant
    @assert sum(abs.(kpb_weight .- bvecs_weight)) < 1e-6

    @debug "k+b k" kpb_k
    @debug "k+b b" kpb_b
    @debug "k+b weights" bvecs_weight

    return BVectors(recip_lattice, kpoints, bvecs, bvecs_weight, kpb_k, kpb_b)
end

@doc raw"""
Generate bvectors for all the kpoints.
"""
function get_bvectors(kpoints::Matrix{T}, recip_lattice::Mat3{T}) where {T<:Real}
    shells = search_shells(kpoints, recip_lattice)

    check_b1(shells)

    bvectors = sort_bvectors(shells)

    # print a blank line to separate the following stdout
    println()

    return bvectors
end
