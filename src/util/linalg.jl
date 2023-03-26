using LinearAlgebra

export get_recip_lattice, get_lattice
export orthonorm_lowdin, eyes_U, rotate_M
export isunitary

"""
    imaglog(z)

Return the imaginary part of the logarithm of `z`.
"""
imaglog(z::T) where {T<:Complex} = atan(imag(z), real(z))

"""
    get_recip_lattice(lattice::Mat3)

Return reciprocal lattice.
"""
get_recip_lattice(lattice::Mat3) = 2π * inv(lattice)'

"""
    get_lattice(recip_lattice::Mat3)

Return lattice.
"""
get_lattice(recip_lattice::Mat3) = inv(recip_lattice / (2π))'

"""
    orthonorm_lowdin(U::Matrix{T})

Lowdin orthonormalize a matrix `U` to be (semi-)unitary.

If `U` is a matrix with orthogonal columns and `V` a non-singular matrix,
then Lowdin-orthogonalizing `U*V` is equivalent to computing `U*orthonorm_lowdin(V)`.
"""
function orthonorm_lowdin(U::Matrix{T}) where {T<:Union{Complex,Real}}
    A, S, B = svd(U)
    # @assert U ≈ A * Diagonal(S) * B'
    return A * B'
end

"""
    orthonorm_lowdin(U::Array{T,3})

Lowdin orthonormalize a series of matrices `M`.
"""
function orthonorm_lowdin(U::Array{T,3}) where {T<:Union{Complex,Real}}
    n_kpts = size(U, 3)

    A = similar(U)

    for ik in 1:n_kpts
        A[:, :, ik] .= orthonorm_lowdin(U[:, :, ik])
    end

    return A
end

function orthonorm_cholesky(U)
    return U / chol(U'U)
end

"""
    fix_global_phase(W::AbstractArray)

Return a factor to fix the global phase of wavefunction,
such that the point having max norm is real.

# Arguments
- `W`: usually `size(W) = nx * ny * nz`

!!! note

    This follows the same logic as `Wannier90` when computing the ratio for real space WFs.
"""
function fix_global_phase(W::AbstractArray)
    f = 1.0 + 0.0im
    # I use abs2 and findmax (returns the 1st maximum)
    # to exactly reproduce W90 behavior
    m, idx = findmax(abs2, W)
    if m > 0
        f = conj(W[idx]) / sqrt(m)
    end
    return f
end

"""
    compute_imre_ratio(W::AbstractArray)

Compute Im/Re ratio of the wavefunction.

!!! note

    This follows the same logic as `Wannier90` when computing the ratio for real space WFs.
"""
function compute_imre_ratio(W::AbstractArray)
    # only calculate real >= 0.01 elements, same as W90
    V = W[abs.(real(W)) .>= 0.01]
    if isempty(V)
        return 0.0
    end
    r = maximum(abs.(imag(V) ./ real(V)))
    return r
end

"""
Power of a unitary (or at least, normal) matrix `U`.
"""
# TODO cleanup this, seems not used anymore
function powm(U::AbstractMatrix{T}, p::F) where {T<:Union{Complex,Real},F<:Real}
    # Workaround, eigen incompatible with lazy adjoint.
    d, V = eigen(Matrix(U))

    V = orthonorm_lowdin(V)
    # accuracy = norm(V * Diagonal(d) * V' - U)
    # @assert accuracy < 1e-10

    return V * Diagonal(d .^ p) * V'
end

"""
    rotate_gauge(O::Array{T,3}, U::Array{T,3})

Rotate the gauge of the operator `O`.

I.e., ``U^{\\dagger} O U``.
"""
function rotate_gauge(O::Array{T,3}, U::Array{T,3}) where {T<:Number}
    n_bands, n_wann, n_kpts = size(U)
    size(O) != (n_bands, n_bands, n_kpts) &&
        error("O must have size (n_bands, n_bands, n_kpts)")

    O1 = similar(O, n_wann, n_wann, n_kpts)

    for ik in 1:n_kpts
        O1[:, :, ik] .= U[:, :, ik]' * O[:, :, ik] * U[:, :, ik]
    end

    return O1
end

"""
    eyes_U(T::Type, n_wann::Int, n_kpts::Int)

Return a series of indentity matrices of type `T` and size `n_wann * n_wann * n_kpts`.
"""
function eyes_U(T::Type, n_wann::Int, n_kpts::Int)
    U = zeros(T, n_wann, n_wann, n_kpts)
    Iₖ = diagm(0 => ones(n_wann))

    for ik in 1:n_kpts
        U[:, :, ik] = Iₖ
    end

    return U
end

"""
    eyes_U(T::Type, n_bands::Int, n_wann::Int, n_kpts::Int)

Return a series of indentity matrices of type `T` and size `n_bands * n_wann * n_kpts`.
"""
function eyes_U(T::Type, n_bands::Int, n_wann::Int, n_kpts::Int)
    U = zeros(T, n_bands, n_wann, n_kpts)
    n = min(n_bands, n_wann)
    Iₖ = diagm(0 => ones(n))

    for ik in 1:n_kpts
        U[1:n, 1:n, ik] = Iₖ
    end

    return U
end

"""
    rotate_U(U::Array{T,3}, V::Array{T,3})

Rotate the gauge matrices `U` by `V`.

I.e., for each kpoint ``\\bm{k}``, ``U_{\\bm{k}} V_{\\bm{k}}``.
"""
function rotate_U(U::Array{T,3}, V::Array{T,3}) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(U)
    size(V)[[1, 3]] != (n_wann, n_kpts) && error("V must be a n_wann * ? * n_kpts matrix")
    m = size(V, 2)

    U1 = similar(U, n_bands, m, n_kpts)

    for ik in 1:n_kpts
        U1[:, :, ik] .= U[:, :, ik] * V[:, :, ik]
    end

    return U1
end

"""
    rotate_M(M::Array{T,4}, kpb_k::Matrix{Int}, U::Array{T,3})

Rotate `mmn` matrices according to gauge `U`.

i.e., for each kpoint ``\\bm{k}``,
``U_{\\bm{k}+\\bm{b}}^{\\dagger} M_{\\bm{k},\\bm{b}} U_{\\bm{k}}``.
"""
@views function rotate_M(
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
    isunitary(U::AbstractArray{T,3}; atol=1e-10)

Check if matrix is unitary or semi-unitary for all the kpoints?

I.e. does it have orthogonal columns?
"""
function isunitary(U::AbstractArray{T,3}; atol::Real=1e-10) where {T<:Number}
    n_bands, n_wann, n_kpts = size(U)

    for ik in 1:n_kpts
        Uₖ = @view U[:, :, ik]
        if norm(Uₖ' * Uₖ - I) > atol
            @debug "not unitary" ik norm(Uₖ' * Uₖ - I)
            return false
        end
    end
    return true
end

"""
    get_projectability(U::AbstractArray{T,3})

Return projectability of each kpoint.
"""
function get_projectability(U::AbstractArray{T,3}) where {T<:Number}
    n_bands, n_wann, n_kpts = size(U)
    P = zeros(real(T), n_bands, n_kpts)
    for ik in 1:n_kpts
        p = U[:, :, ik] * U[:, :, ik]'
        P[:, ik] = real(diag(p))
    end
    return P
end

"""
    findvector(predicate::Function, v::AbstractVector, M::AbstractMatrix)

Find index of vector in the columns of a matrix.

# Arguments
- `predicate`: comparison function
- `v`: the vector to be found
- `M`: the matrix to be searched, its column will be compared to `v` by `predicate` function
"""
function findvector(predicate::Function, v::AbstractVector, M::AbstractMatrix)
    for (i, col) in enumerate(eachcol(M))
        predicate(v, col) && return i
    end
    error("$v not found in array!")
    return nothing
end

"""
    rand_unitary(T::Type, m::Int, n::Int)

Generate a random (semi-)unitary matrix using Lowdin orthonormalization.

# Arguments
- `T`: the type of the matrix, e.g., `ComplexF64`, `Float64`
- `m`: number of rows
- `n`: number of columns
"""
function rand_unitary(T::Type, m::Int, n::Int)
    M = randn(T, m, n)
    N = orthonorm_lowdin(M)
    return N
end

"""
    rand_unitary(T::Type, m::Int)

Generate a random unitary matrix using Lowdin orthonormalization.

# Arguments
- `T`: the type of the matrix
- `m`: number of rows (= number of columns)
"""
rand_unitary(T::Type, m::Int) = rand_unitary(T, m, m)

"""
    rand_unitary(T::Type, m::Int, n::Int, k::Int)

Generate a series of random (semi-)unitary matrix using Lowdin orthonormalization.

The returned `M[:, :, ik]` is (semi-)unitary for all `ik = 1:k`.

# Arguments
- `T`: the type of the matrix
- `m`: number of rows
- `n`: number of columns
- `k`: number of matrices
"""
function rand_unitary(T::Type, m::Int, n::Int, k::Int)
    M = zeros(T, m, n, k)
    for ik in 1:k
        M[:, :, ik] = rand_unitary(T, m, n)
    end
    return M
end
