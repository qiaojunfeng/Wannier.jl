import LinearAlgebra as LA


imaglog(z::T) where {T<:Complex} = atan(imag(z), real(z))


function get_recip_lattice(lattice::Mat3)
    2pi * inv(lattice)'
end


@doc raw"""
Computes overlap between two neighboring kpoints

size(M) = (n_bands, n_bands, n_bvecs, n_kpts)
size(kpb_k) = (n_bvecs, n_kpts)
"""
function overlap(
    M::Array{Complex{FT},4},
    kpb_k::Matrix{Int},
    k1::Int,
    k2::Int,
) where {FT<:Real}
    n_bands = size(M, 1)
    n_bvecs = size(M, 3)

    if (k1 == k2)
        return Matrix{Complex{FT}}((1.0 + 0.0im)LA.I, n_bands, n_bands)
    end

    for ib = 1:n_bvecs
        if kpb_k[ib, k1] == k2
            return view(M, :, :, ib, k1)
        end
    end

    error("No neighbors found, k1 = $(k1), k2 = $(k2)")

    # Matrix{Complex{FT}}((1.0 + 0.0im)LA.I, n_bands, n_bands)
end

@doc raw"""
Computes overlap between two neighboring kpoints, and rotate by gauge from A matrices.

size(M) = (n_bands, n_bands, n_bvecs, n_kpts)
size(kpb_k) = (n_bvecs, n_kpts)
size(A) = (n_bands, n_wann, n_kpts)
"""
function overlap(
    M::Array{Complex{FT},4},
    kpb_k::Matrix{Int},
    k1::Int,
    k2::Int,
    A::Array{Complex{FT},3},
) where {FT<:Real}
    A[:, :, k1]' * overlap(M, kpb_k, k1, k2) * A[:, :, k2]
end


"""
Normalize a matrix A to be (semi-)unitary. 
If X is a matrix with orthogonal columns and A a non-singular matrix,
then Lowdin-orthogonalizing X*A is equivalent to computing X*normalize_matrix(A)
"""
function orthonormalize_lowdin(A)
    U, S, V = LA.svd(A)
    @assert A ≈ U * LA.Diagonal(S) * V'
    return U * V'
end


function orthonormalize_cholesky(A)
    ovl = A'A
    return A / LA.chol(ovl)
end
