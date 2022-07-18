using LinearAlgebra

imaglog(z::T) where {T<:Complex} = atan(imag(z), real(z))

get_recip_lattice(lattice::Mat3) = 2π * inv(lattice)'
get_lattice(recip_lattice::Mat3) = inv(recip_lattice / (2π))'

@doc raw"""
Computes overlap between two neighboring kpoints

size(M) = (n_bands, n_bands, n_bvecs, n_kpts)
size(kpb_k) = (n_bvecs, n_kpts)
"""
function overlap(
    M::Array{Complex{FT},4}, kpb_k::Matrix{Int}, k1::Int, k2::Int
) where {FT<:Real}
    n_bands = size(M, 1)
    n_bvecs = size(M, 3)

    if (k1 == k2)
        return Matrix{Complex{FT}}((1.0 + 0.0im)I, n_bands, n_bands)
    end

    for ib in 1:n_bvecs
        if kpb_k[ib, k1] == k2
            return view(M, :, :, ib, k1)
        end
    end

    return error("No neighbors found, k1 = $(k1), k2 = $(k2)")

    # Matrix{Complex{FT}}((1.0 + 0.0im)I, n_bands, n_bands)
end

@doc raw"""
Computes overlap between two neighboring kpoints, and rotate by gauge from A matrices.

size(M) = (n_bands, n_bands, n_bvecs, n_kpts)
size(kpb_k) = (n_bvecs, n_kpts)
size(A) = (n_bands, n_wann, n_kpts)
"""
function overlap(
    M::Array{Complex{FT},4}, kpb_k::Matrix{Int}, k1::Int, k2::Int, A::Array{Complex{FT},3}
) where {FT<:Real}
    return A[:, :, k1]' * overlap(M, kpb_k, k1, k2) * A[:, :, k2]
end

"""
Normalize a matrix A to be (semi-)unitary.
If X is a matrix with orthogonal columns and A a non-singular matrix,
then Lowdin-orthogonalizing X*A is equivalent to computing X*normalize_matrix(A)
"""
function orthonorm_lowdin(A::Matrix{T}) where {T<:Union{Complex,Real}}
    U, S, V = svd(A)
    # @assert A ≈ U * Diagonal(S) * V'
    return U * V'
end

function orthonorm_lowdin(A::Array{T,3}) where {T<:Union{Complex,Real}}
    n_kpts = size(A, 3)

    B = similar(A)

    for ik in 1:n_kpts
        B[:, :, ik] .= orthonorm_lowdin(A[:, :, ik])
    end

    return B
end

function orthonorm_cholesky(A)
    return A / chol(A'A)
end

function fix_global_phase!(A::Array{T,3}) where {T<:Complex}
    n_wann = size(A, 2)

    for i in 1:n_wann
        imax = indmax(abs.(A[:, i, 1]))
        a = A[imax, i, 1]
        @assert abs(a) > 1e-2
        A[:, i, :] *= conj(a / abs(a))
    end

    return nothing
end

"""
Power of a unitary (or at least, normal) matrix A
"""
function powm(A::AbstractMatrix{T}, p::F) where {T<:Union{Complex,Real},F<:Real}
    # Workaround, eigen incompatible with lazy adjoint.
    d, V = eigen(Matrix(A))

    V = orthonorm_lowdin(V)
    # accuracy = norm(V * Diagonal(d) * V' - A)
    # @assert accuracy < 1e-10

    return V * Diagonal(d .^ p) * V'
end

function rotate_gauge(O::Array{T,3}, A::Array{T,3}) where {T<:Number}
    n_bands, n_wann, n_kpts = size(A)
    size(O) != (n_bands, n_bands, n_kpts) &&
        error("O must have size (n_bands, n_bands, n_kpts)")

    O1 = similar(O, n_wann, n_wann, n_kpts)

    for ik in 1:n_kpts
        O1[:, :, ik] .= A[:, :, ik]' * O[:, :, ik] * A[:, :, ik]
    end

    return O1
end

function eyes_amn(T::Type, n_wann::Int, n_kpts::Int)
    A = zeros(T, n_wann, n_wann, n_kpts)
    Iₖ = diagm(0 => ones(n_wann))

    for ik in 1:n_kpts
        A[:, :, ik] = Iₖ
    end

    return A
end

function eyes_amn(T::Type, n_bands::Int, n_wann::Int, n_kpts::Int)
    A = zeros(T, n_bands, n_wann, n_kpts)
    n = min(n_bands, n_wann)
    Iₖ = diagm(0 => ones(n))

    for ik in 1:n_kpts
        A[1:n, 1:n, ik] = Iₖ
    end

    return A
end

"""
Rotate MMN matrices according to gauge U.
"""
@views function rotate_mmn(
    M::Array{T,4}, kpb_k::Matrix{Int}, U::Array{T,3}
) where {T<:Complex}
    n_bands, n_wann = size(U)
    n_kpts = size(M, 4)
    n_bvecs = size(M, 3)

    n_bands != size(M, 1) && error("incompatible n_bands")

    # Fill MMN
    N = similar(M, n_wann, n_wann, n_bvecs, n_kpts)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ik2 = kpb_k[ib, ik]

            U₁ = U[:, :, ik]
            U₂ = U[:, :, ik2]

            N[:, :, ib, ik] = U₁' * M[:, :, ib, ik] * U₂
        end
    end

    return N
end

"""
Is matrix unitary or semi-unitary for all the kpoints?
i.e. does it have orthogonal columns?
"""
function isunitary(A::AbstractArray{T,3}; atol::Real=1e-10) where {T<:Number}
    n_bands, n_wann, n_kpts = size(A)

    for ik in 1:n_kpts
        Aₖ = @view A[:, :, ik]
        if norm(Aₖ' * Aₖ - I) > atol
            @debug "not unitary" ik norm(Aₖ' * Aₖ - I)
            return false
        end
    end
    return true
end
