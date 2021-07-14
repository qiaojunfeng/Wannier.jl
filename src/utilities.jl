module Utilities

import LinearAlgebra as LA

function get_recipcell(unit_cell)
    return 2pi * inv(unit_cell)'
end

# Normalize a matrix A to be (semi-)unitary. 
# If X is a matrix with orthogonal columns and A a non-singular matrix, then Lowdin-orthogonalizing X*A is equivalent to computing X*normalize_matrix(A)
function orthonormalize_lowdin(A)
    U, S, V = LA.svd(A)
    @assert A ≈ U * LA.Diagonal(S) * V'
    return U * V'
end

function orthonormalize_cholesky(A)
    ovl = A'A
    return A / LA.chol(ovl)
end

# Computes overlap between two neighboring kpoints
function overlap(data, k1, k2)
    # fix useful if N3 = 1 (e.g. for 2D models)
    if (k1 == k2)
        return Matrix((1.0 + 0.0im)I, data.num_bands, data.num_bands)
    end
    for ib = 1:data.num_bvecs
        if data.kpbs[ib, k1] == k2
            return view(data.mmn, :, :, ib, k1)
        end
    end
    error("No neighbors found, k1 = $(k1), k2 = $(k2)")
    return Matrix((1.0 + 0.0im)I, data.num_bands, data.num_bands)
end

function get_projectability(A)
    num_bands, num_wann, num_kpts = size(A)
    proj = zeros(num_bands, num_kpts)
    for ik = 1:num_kpts
        p = A[:,:,ik] * A[:,:,ik]'
        proj[:,ik] = real(LA.diag(p))
    end
    return proj
end

end
