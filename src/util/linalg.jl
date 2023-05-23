using LinearAlgebra

export orthonorm_lowdin, eyes_U, rotate_M
export isunitary

"""
    imaglog(z)

Return the imaginary part of the logarithm of `z`.
"""
imaglog(z::T) where {T<:Complex} = atan(imag(z), real(z))

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
orthonorm_lowdin(U::Vector{Matrix{T}}) where {T<:Union{Complex,Real}} = orthonorm_lowdin.(U)

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
function rotate_gauge(O::Vector{Matrix{T}}, U::Vector{Matrix{T}}) where {T<:Number}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    size(O[1], 1), size(O[1], 2), length(O) != (n_bands, n_bands, n_kpts) &&
        error("O must have size (n_bands, n_bands, n_kpts)")

    return map(zip(O, U)) do (o, u)
        u' * o * u
    end
end

"""
    eyes_U(T::Type, n_wann::Int, n_kpts::Int)

Return a series of indentity matrices of type `T` and size `n_wann * n_wann * n_kpts`.
"""
eyes_U(::Type{T}, n_wann::Int, n_kpts::Int) where T = 
    [diagm(0 => ones(T, n_wann)) for i = 1:n_kpts]

"""
    eyes_U(T::Type, n_bands::Int, n_wann::Int, n_kpts::Int)

Return a series of indentity matrices of type `T` and size `n_bands * n_wann * n_kpts`.
"""
function eyes_U(T::Type, n_bands::Int, n_wann::Int, n_kpts::Int)
    map(1:n_kpts) do
        U = zeros(T, n_bands, n_wann)
        n = min(n_bands, n_wann)
        U[1:n, 1:n] .= 1
        return U
    end
end

@doc raw"""
    rotate_U(U, V)

Rotate the gauge matrices `U` by `V`.

For each kpoint ``\bm{k}``, return ``U_{\bm{k}} V_{\bm{k}}``.

# Arguments
- `U`: a series of gauge matrices, usually `size(U) = n_bands * n_wann * n_kpts`
- `V`: a series of gauge matrices, usually `size(V) = n_wann * n_wann * n_kpts`
"""
function rotate_U(U::AbstractVector, V::AbstractVector)
    n_kpts = length(U)
    n_bands, n_wann = size(U[1])
    
    return map(zip(U, V)) do (u, v)
        u * v
    end
end

"""
    rotate_M(M::Array{T,4}, kpb_k::Matrix{Int}, U::Array{T,3})

Rotate `mmn` matrices according to gauge `U`.

i.e., for each kpoint ``\\bm{k}``,
``U_{\\bm{k}+\\bm{b}}^{\\dagger} M_{\\bm{k},\\bm{b}} U_{\\bm{k}}``.
"""
@views function rotate_M(
    M::Vector{Array{T,3}}, kpb_k, U::Vector{Matrix{T}}
) where {T<:Complex}

    n_bands, n_wann = size(U[1])
    n_kpts = length(M)
    n_bvecs = size(M[1], 3)

    n_bands != size(M[1], 1) && error("incompatible n_bands")

    # Fill MMN
    N = [similar(M[1], n_wann, n_wann, n_bvecs) for i = 1:n_kpts]

    @views @inbounds for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ik2 = kpb_k[ik][ib]
            U₁ = U[ik]
            U₂ = U[ik2]

            N[ik][:, :, ib] .= U₁' * M[ik][:, :, ib] * U₂
        end
    end

    return N
end

"""
    isunitary(U::AbstractVector{AbstractMatrix{T}}; atol=1e-10)

Check if matrix is unitary or semi-unitary for all the kpoints?

I.e. does it have orthogonal columns?
"""
function isunitary(U::AbstractVector{AbstractMatrix{T}}; atol::Real=1e-10) where {T<:Number}
    n_bands, n_wann, n_kpts = size(U)

    for ik in 1:n_kpts
        Uₖ = U[ik]
        if norm(Uₖ' * Uₖ - I) > atol
            @debug "not unitary" ik norm(Uₖ' * Uₖ - I)
            return false
        end
    end
    return true
end

"""
    get_projectability(U::AbstractVector{Abs = size(U[1])
    actMatrix{T}})

Return projectability of each kpoint.
"""
function get_projectability(U::AbstractVector{AbstractMatrix{T}}) where {T<:Number}
    map(U) do u
        p = u * u'
        return real(diag(p))
    end
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
rand_unitary(T::Type, m::Int, n::Int, k::Int) = [rand_unitary(T, m, n) for i = 1:k]
