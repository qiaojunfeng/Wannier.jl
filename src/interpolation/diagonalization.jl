using FastLapackInterface: HermitianEigenWs, decompose!

struct _HermitianDiagonalizationWorkspace{W}
    lapack::W
end

function _HermitianDiagonalizationWorkspace(
        ::Type{T}, number_wannier::Integer
    ) where {T <: Complex}
    matrix = zeros(T, number_wannier, number_wannier)
    return _HermitianDiagonalizationWorkspace(HermitianEigenWs(matrix))
end

function _diagonalize_hermitian_batch!(
        eigenvalues::AbstractMatrix,
        eigenvectors::AbstractArray{<:Complex, 3},
        matrices::AbstractArray{<:Complex, 3},
        workspace::_HermitianDiagonalizationWorkspace,
    )
    size(eigenvectors) == size(matrices) ||
        throw(DimensionMismatch("matrix and eigenvector batches differ"))
    size(eigenvalues) == (size(matrices, 1), size(matrices, 3)) ||
        throw(DimensionMismatch("matrix and eigenvalue batches differ"))

    for batch_index in axes(matrices, 3)
        vectors = view(eigenvectors, :, :, batch_index)
        copyto!(vectors, view(matrices, :, :, batch_index))
        decompose!(workspace.lapack, 'V', 'A', 'U', vectors, 0.0, 0.0, 0, 0, 1.0e-16)
        copyto!(view(eigenvalues, :, batch_index), workspace.lapack.w)
        copyto!(vectors, workspace.lapack.Z)
    end
    return eigenvalues, eigenvectors
end
