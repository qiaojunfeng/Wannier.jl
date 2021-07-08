module BVectors

import LinearAlgebra as LA
import NearestNeighbors as NN

mutable struct BVectorShells
    # reciprocal cell, 3 * 3
    recip_cell::Array{Float64,2}

    # kpoints array, reduced coordinates, 3 * num_kpts
    kpts::Array{Float64,2}

    # number of shells
    num_shells::Int

    # bvectors of each shell, cartesian coordinates, 3 * max_multiplicity * num_shells
    vecs::Array{Float64,3}

    # multiplicity of each shell, num_shells
    multis::Array{Int,1}

    # weight of each shell, num_shells
    weights::Array{Float64,1}
end

"""
kpts and returned supercell are in reduced coord
"""
function generate_supercell(kpts)
    # The total space to be searched
    supercell_size = 5

    num_cell = (2 * supercell_size + 1)^3
    num_kpts = size(kpts, 2)
    supercell = zeros(3, num_cell * num_kpts)
    supercell_idx = zeros(3, num_cell * num_kpts)
    counter = 1
    # Note: the order of which index increases the fastest is crucial,
    # it will determine the order of bvecotors, 
    for i = -supercell_size:supercell_size
        for j = -supercell_size:supercell_size
            for k = -supercell_size:supercell_size
                supercell_idx[:, counter:(counter + num_kpts - 1)] .= [i, j, k]
                supercell[:, counter:(counter + num_kpts - 1)] = kpts .+ [i, j, k]
                counter += num_kpts
            end
        end
    end

    return supercell, supercell_idx
end

"""
check if the columns of matrix m and columns of matrix n are parallel
"""
function check_parallel(m, n)
    @assert size(m, 1) == size(n, 1) == 3

    eps = 1e-5
    
    result = fill(false, size(m, 2), size(n, 2))
    for (i, col1) in enumerate(eachcol(m))
        for (j, col2) in enumerate(eachcol(n))
            p = LA.cross(col1, col2)
            if all(isapprox.(0, p; atol=eps))
                result[i, j] = true
            end
        end
    end
    return result
end

"""
try to guess weights from MV1997 Eq. B1

    bvecs: 3 * max_multi * num_shell
    multis: num_shell
    return: same size as input, but only selected shells which satisfy B1 condition
"""
function calculate_b1_weights(bvecs, multis)
    eps = 1e-8

    num_shell = length(multis)

    # only compare the upper triangular part of bvec * bvec', 6 elements
    bmat = zeros(6, num_shell)
    # return the upper triangular part of a matrix as a vector
    triu2vec(m) = m[LA.triu!(trues(size(m)), 0)]
    W = zeros(num_shell)
    ishell = 1
    # triu2vec(I) = [1 0 1 0 0 1]
    triu_I = triu2vec(LA.diagm([1,1,1]))
    while ishell <= num_shell
        bmat[:, ishell] = triu2vec(bvecs[:,1:multis[ishell],ishell] * bvecs[:,1:multis[ishell],ishell]')
        # Solve equation bmat * W = triu_I
        # size(bmat) = (6, ishell), W is diagonal matrix of size ishell
        # bmat = U * S * V' -> W = V * S^-1 * U' * triu2vec(I)
        U, S, V = LA.svd(bmat[:, 1:ishell])
        # @debug "S" ishell bmat[:,ishell]' S
        if all(S .> eps)
            W[1:ishell] = V * LA.Diagonal(S)^-1 * U' * triu_I
            if isapprox(bmat[:, 1:ishell] * W[1:ishell], triu_I; atol=eps)
                break
            end
        end
        ishell += 1
    end
    @assert ishell != num_shell + 1 "Did not find enough shells!"
    num_shell = ishell

    # resize
    multis = multis[1:num_shell]
    max_multi = maximum(multis)
    weights = W[1:num_shell]
    bvecs = bvecs[:, 1:max_multi, 1:num_shell]

    return bvecs, multis, weights
end

"""
kpts: reduced coordinates
"""
function search_shells(kpts, recip_cell)
    # Number of nearest-neighbors to be returned
    search_neighbors = 100
    # precision to select a shell (equal distance to the point)
    eps = 1e-5
    # Max number of stencils in one shell
    max_multiplicities = 40
    # Max number of nearest-neighbor shells
    max_shells = round(Int, search_neighbors / max_multiplicities) # 5

    supercell, _ = generate_supercell(kpts)

    # To cartesian coordinates
    supercell_cart = recip_cell * supercell
    kpt_cart = recip_cell * kpts[:, 1] # use the 1st kpt to search bvectors

    kdtree = NN.KDTree(supercell_cart)
    idxs, dists = NN.knn(kdtree, kpt_cart, search_neighbors, true)

    # activate debug info with: JULIA_DEBUG=Main julia
    # @debug "KDTree nearest neighbors" dists
    # @debug "KDTree nearest neighbors" idxs

    multis = zeros(Int, max_shells) # multiplicities per shell
    nearest_kpts = zeros(Int, max_multiplicities, max_shells)

    # The 1st result is the search point itself, dist = 0
    ineigh = 2
    ishell = 1
    while (ineigh <= search_neighbors) && (ishell <= max_shells)
        # use the 1st kpoint to find bvector shells & weights
        eqdist_idxs = findall(x -> isapprox(x, dists[ineigh]; atol=eps), dists)
        multis[ishell] = length(eqdist_idxs)
        if multis[ishell] >= max_multiplicities
            # skip large multiplicity shells
            break
        end
        nearest_kpts[1:multis[ishell], ishell] = idxs[eqdist_idxs]
        ineigh += multis[ishell]
        ishell += 1
    end


    # resize
    num_shell = ishell - 1
    multis = multis[1:num_shell]
    max_multi = maximum(multis)
    bvecs = zeros(3, max_multi, num_shell) # in Cartesian coord
    for ishell = 1:num_shell
        kpb_cart = supercell_cart[:, nearest_kpts[1:multis[ishell], ishell]]
        bvecs[:, 1:multis[ishell], ishell] = kpb_cart .- kpt_cart
    end
    @debug "Found bvector shells" multis
    @debug "Found bvector shells" bvecs
    
    # remove shells with parallel bvectors
    keepshells = collect(1:num_shell)
    for ishell = 2:num_shell
        hasparallel = false
        for jshell = 1:(ishell - 1)
            if !(jshell in keepshells)
                continue
            end
            if hasparallel
                break
            end
            p = check_parallel(bvecs[:, 1:multis[jshell], jshell], bvecs[:, 1:multis[ishell], ishell])
            # @debug "check_parallel($jshell, $ishell)" p
            if any(p)
                hasparallel = true
                break
            end
        end
        if hasparallel
            splice!(keepshells, ishell)
        end
    end
    num_shell = length(keepshells)
    multis = multis[keepshells]
    max_multi = maximum(multis)
    bvecs = bvecs[:, 1:max_multi, keepshells]

    # @debug "After check_parallel bvector shells" multis
    # @debug "After check_parallel bvector shells" bvecs

    # sort bvectors

    bvecs, multis, weights = calculate_b1_weights(bvecs, multis)

    @info "BVector shell vectors" bvecs
    @info "BVector shell multiplicities" multis
    @info "BVector shell weights" weights

    return BVectorShells(recip_cell, kpts, num_shell, bvecs, multis, weights)
end

function check_b1(shells)
    # Check B1
    mat = zeros(3, 3)
    for ish = 1:length(shells.multis)
        bvec_ish = shells.vecs[:, 1:shells.multis[ish], ish]
        mat += shells.weights[ish] * bvec_ish * bvec_ish'
    end
    @debug "Bvector sum" mat
    @assert isapprox(mat, LA.I)
    println("Finite difference condition satisfied")
end

"""
flatten shell vectors
return: 
bvecs: 3 * num_bvecs
bvecs_weight: num_bvecs
"""
function flatten_shells(shells)
    num_bvecs = sum(shells.multis)
    bvecs = zeros(3, num_bvecs)
    bvecs_weight = zeros(num_bvecs)
    counter = 1
    for ishell = 1:length(shells.multis)
        imulti = shells.multis[ishell]
        bvecs[:, counter:(counter + imulti - 1)] = shells.vecs[:, 1:imulti, ishell]
        bvecs_weight[counter:(counter + imulti - 1)] .= shells.weights[ishell]
        counter += imulti
    end
    return bvecs, bvecs_weight
end

# Generate bvectors for all the kpoints & sort (to be consistent with wannier90)
function generate_bvectors(kpts, recip_cell)

    shells = search_shells(kpts, recip_cell)

    check_b1(shells)

    num_kpts = size(kpts, 2)

    kpts_cart = recip_cell * shells.kpts
    # Calculate distances of supercells to original cell
    _, supercell_idx = generate_supercell(zeros(3))

    # Firstly, sort by length of bvectors: nearest k+b goes first
    # Secondly, sort by supercell: k+b which is in the same cell as k goes first
    # Thirdly, sort by the index (in the supercell array) of supercell: the larger-index one goes first
    #     Note this also relies on the order of supercell_cart, see generate_supercell()
    # Fourthly, sort by k+b index: the smaller-index k+b point goes first
    inv_recip = inv(recip_cell)
    
    function isequiv(v1, v2)
        eps = 1e-5
        d = v1 - v2
        d .-= round.(d)
        return all(isapprox.(d, 0; atol=eps))
    end

    function findvector(predicate::Function, v, vecs)
        for (i, col) in enumerate(eachcol(vecs))
            predicate(v, col) && return i
        end
        throw(ErrorException("$v not found in array!"))
    end

    """
    kb1, kb2, k: cartesian coordinates
    """
    function bvec_isless(kb1, kb2, k)
        # Float64 comparsion is tricky, esp. for norm
        eps = 1e-5
        
        if isapprox(kb1, kb2; atol=eps)
            return true
        end

        normkb1 = LA.norm(kb1 - k)
        normkb2 = LA.norm(kb2 - k)
        if abs(normkb1 - normkb2) < eps
            # using outter scope variable: inv_recip, kpts, kpts_cart
            # bring back to reduced coord, easier to compare
            kb1_reduced = inv_recip * kb1
            kb2_reduced = inv_recip * kb2
            kb1_equiv_idx = findvector(isequiv, kb1_reduced, kpts)
            kb2_equiv_idx = findvector(isequiv, kb2_reduced, kpts)
            celldisp_kb1 = kb1 - kpts_cart[:, kb1_equiv_idx]
            celldisp_kb2 = kb2 - kpts_cart[:, kb2_equiv_idx]
            # @debug "jq" celldisp_kb1, celldisp_kb2
            normkb1 = LA.norm(celldisp_kb1)
            normkb2 = LA.norm(celldisp_kb2)
            # @debug "jq" normkb1, normkb2
            if abs(normkb1 - normkb2) < eps
                # using outter scope variable: supercell_idx
                celldisp_kb1_reduced = round.(Int, inv_recip * celldisp_kb1)
                celldisp_kb2_reduced = round.(Int, inv_recip * celldisp_kb2)
                idx1 = findvector(==, celldisp_kb1_reduced, supercell_idx)
                idx2 = findvector(==, celldisp_kb2_reduced, supercell_idx)
                # @debug "jq", idx1, idx2
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
                        throw(ErrorException("you are comparing equivalent k+b points, this should not happen!"))
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

    # FIX: the order of b vectors is not the same as wannier90
    # k = kpts_cart[:,1]
    # kb1 = recip_cell * (kpts[:,9] +[0 0 0]')
    # kb2 = recip_cell * (kpts[:,65] +[0 0 0]')
    # @debug "islesss" bvec_isless(kb1, kb2, k)
    # throw(ErrorException)

    bvecs, bvecs_weight = flatten_shells(shells)
    num_bvecs = length(bvecs_weight)

    # find k+b indexes
    kpbs = zeros(Int, num_bvecs, num_kpts)
    kpbs_weight = zeros(num_bvecs, num_kpts)
    kpbs_disp = zeros(Int, 3, num_bvecs, num_kpts)
    sorted_idx = zeros(Int, num_bvecs)
    # I have to construct a Vector of Vector such that sortperm could work
    ikpb = Vector{Vector{Float64}}(undef, num_bvecs)
    for ik = 1:num_kpts
        k_cart = kpts_cart[:, ik]
        for ib = 1:num_bvecs
            ikpb[ib] = k_cart .+ bvecs[:, ib]
        end
        sorted_idx[:] = sortperm(ikpb; lt=(x, y) -> bvec_isless(x, y, k_cart))
        for ib = 1:num_bvecs
            ikpb_reduced = inv_recip * ikpb[sorted_idx[ib]]
            ikpb_equiv_idx = findvector(isequiv, ikpb_reduced, kpts)
            kpbs[ib, ik] = ikpb_equiv_idx
            celldisp = ikpb_reduced - kpts[:, ikpb_equiv_idx]
            kpbs_disp[:, ib, ik] = round.(Int, celldisp)
            kpbs_weight[ib, ik] = bvecs_weight[sorted_idx[ib]]
        end
    end
    @debug "k+b" kpbs
    @debug "k+b displacements" kpbs_disp
    @debug "k+b weights" kpbs_weight

    return kpbs, kpbs_disp, kpbs_weight
end

"""
print Wannier90 nnkp file bvectors in cartesian coordinates
"""
function print_w90_nnkp(w90nnkp)
    for k = 1:w90nnkp.num_kpts
        kpt = w90nnkp.kpoints[:, k]
        println("k = $kpt")
        for b = 1:w90nnkp.num_bvecs
            nnkp = w90nnkp.nnkpts[:, b, k]
            kpb = w90nnkp.kpoints[:, nnkp[1]] + nnkp[2:end]
            bvec = kpb - kpt
            bvec_cart = w90nnkp.recip_lattice * bvec
            println(bvec_cart)
        end
    end
end

end
