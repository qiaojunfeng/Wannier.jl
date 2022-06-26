
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
