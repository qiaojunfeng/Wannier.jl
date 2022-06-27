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
function orthonorm_lowdin(A::Matrix{T}) where {T<:Union{Complex,Real}}
    U, S, V = LA.svd(A)
    # @assert A ≈ U * LA.Diagonal(S) * V'
    U * V'
end


function orthonorm_lowdin(A::Array{T,3}) where {T<:Union{Complex,Real}}
    n_kpts = size(A, 3)

    B = similar(A)

    for ik = 1:n_kpts
        B[:, :, ik] .= orthonorm_lowdin(A[:, :, ik])
    end

    B
end


function orthonorm_cholesky(A)
    A / LA.chol(A'A)
end


function fix_global_phase!(A::Array{T,3}) where {T<:Complex}
    n_wann = size(A, 2)

    for i = 1:n_wann
        imax = indmax(abs.(A[:, i, 1]))
        a = A[imax, i, 1]
        @assert abs(a) > 1e-2
        A[:, i, :] *= conj(a / abs(a))
    end

    nothing
end
