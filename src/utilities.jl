module Utilities

using LinearAlgebra

function get_recipcell(unit_cell)
    return 2pi * inv(unit_cell)'
end

# Normalize a matrix A to be (semi-)unitary. 
# If X is a matrix with orthogonal columns and A a non-singular matrix, then Lowdin-orthogonalizing X*A is equivalent to computing X*normalize_matrix(A)
function orthonormalize_lowdin(A)
    U, S, V = svd(A)
    @assert A ≈ U * Diagonal(S) * V'
    return U * V'
end

function orthonormalize_cholesky(A)
    ovl = A'A
    return A / chol(ovl)
end

# Computes overlap between two neighboring kpoints
function overlap(params, k1, k2)
    # fix useful if N3 = 1 (e.g. for 2D models)
    if (k1 == k2)
        return Matrix((1.0 + 0.0im)I, params.num_bands, params.num_bands)
    end
    for ib = 1:params.num_bvecs
        if params.kpbs[ib, k1] == k2
            return view(params.mmn, :, :, ib, k1)
        end
    end
    error("No neighbors found, k1 = $(k1), k2 = $(k2)")
    return Matrix((1.0 + 0.0im)I, params.num_bands, params.num_bands)
end

end
