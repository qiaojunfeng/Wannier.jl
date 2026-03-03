export orthonorm_lowdin,
    identity_gauge,
    zeros_gauge,
    rand_gauge,
    merge_gauge,
    zeros_overlap,
    zeros_eigenvalues,
    isunitary,
    isequiv

"""
    $(SIGNATURES)

Return the imaginary part of the logarithm of `z`.
"""
imaglog(z::Complex) = atan(imag(z), real(z))

"""
    $(SIGNATURES)

Lowdin orthonormalize a matrix `U` to be (semi-)unitary.

If `U` is a matrix with orthogonal columns and `V` a non-singular matrix,
then Lowdin-orthogonalizing `U*V` is equivalent to computing `U*orthonorm_lowdin(V)`.
"""
function orthonorm_lowdin(U::AbstractMatrix)
    A, S, B = svd(U)
    # @assert U ≈ A * Diagonal(S) * B'
    return A * B'
end

"""
    $(SIGNATURES)

Lowdin orthonormalize a series of matrices `U`.
"""
orthonorm_lowdin(U::AbstractVector) = orthonorm_lowdin.(U)

function orthonorm_cholesky(U)
    return U / chol(U'U)
end

"""
    $(SIGNATURES)

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
    $(SIGNATURES)

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
    $(SIGNATURES)

Allocate gauge matrix `U` filled with identity matrices.

The `U` can be accessed by `U[ik][m, n]`, where `ik`, `m`, `n` are
the indices of kpoints, bands, and WFs, respectively.

# Arguments
- `T`: the type of the matrix elements
- `nkpts`: number of kpoints
- `nwann`: number of Wannier functions
"""
function identity_gauge(T::Type, nkpts::Integer, nwann::Integer)
    return map(1:nkpts) do _
        diagm(0 => ones(T, nwann))
    end
end

"""
    zeros_gauge(T, nkpts, nbands, nwann)
    zeros_gauge(T, nkpts, nwann)

Allocate gauge matrix `U` filled with zeros.

See also [`identity_gauge`](@ref).
"""
function zeros_gauge end

function zeros_gauge(T::Type, nkpts::Integer, nbands::Integer, nwann::Integer)
    return map(1:nkpts) do _
        zeros(T, nbands, nwann)
    end
end

function zeros_gauge(T::Type, nkpts::Integer, nwann::Integer)
    return zeros_gauge(T, nkpts, nwann, nwann)
end

"""
    $(SIGNATURES)

Merge the two sets of gauge matrices `U` and `V`.

For each kpoint ``\\mathbf{k}``, return ``U_{\\mathbf{k}} V_{\\mathbf{k}}``.

# Arguments
- `U`: a series of gauge matrices
- `V`: a series of gauge matrices
"""
function merge_gauge(U::AbstractVector, V::AbstractVector)
    return map(zip(U, V)) do (u, v)
        u * v
    end
end

"""
    rand_gauge(T, nwann)
    rand_gauge(T, nbands, nwann)
    rand_gauge(T, nkpts, nbands, nwann)

Generate random (semi-)unitary matrices using Lowdin orthonormalization.

The returned `U` can be accessed by `U[ik][m, n]`, where `ik`, `m`, `n` are
the indices of kpoints, bands, and WFs, respectively, and each `U[ik]` is
(semi-)unitary.

# Arguments
- `T`: the type of the matrix elements, e.g., `ComplexF64`, `Float64`
- `nkpts`: number of kpoints
- `nbands`: number of bands
- `nwann`: number of Wannier functions
"""
function rand_gauge end

@inline function rand_gauge(T::Type, nbands::Integer, nwann::Integer)
    return orthonorm_lowdin(randn(T, nbands, nwann))
end

@inline rand_gauge(T::Type, nwann::Integer) = rand_gauge(T, nwann, nwann)

function rand_gauge(T::Type, nkpts::Integer, nbands::Integer, nwann::Integer)
    return map(1:nkpts) do _
        rand_gauge(T, nbands, nwann)
    end
end

"""
    $(SIGNATURES)

Allocate overlap `M` matrices filled with zeros.

The `M` can be accessed by `M[ik][ib][m, n]`, where `ik`, `ib`, `m`, `n` are
the indices of kpoints, b-vectors, bands, and WFs, respectively.
"""
function zeros_overlap(T::Type, nkpts::Integer, nbvecs::Integer, nbands::Integer)
    return map(1:nkpts) do _
        map(1:nbvecs) do _
            zeros(T, nbands, nbands)
        end
    end
end

"""
    $(SIGNATURES)

Allocate eigenvalues matrices filled with zeros.

The returned `E` can be accessed by `E[ik][m]`, where `ik`, `m` are the indices
of kpoints and bands, respectively.
"""
function zeros_eigenvalues(T::Type, nkpts::Integer, nwann::Integer)
    return [zeros(T, nwann) for _ in 1:nkpts]
end

"""
    $(SIGNATURES)

Check if matrix is unitary or semi-unitary for all the kpoints.

I.e. if it has orthogonal columns.
"""
function isunitary(U::AbstractVector; atol=1e-10)
    map(U) do u
        if norm(u' * u - I) > atol
            @debug "not unitary" ik norm(u' * u - I)
            return false
        end
    end
    return true
end

"""
    $(SIGNATURES)

Return projectability of each kpoint.

The projections (gauge matrices) are <ψ|ϕ> where |ψ> is the Bloch state and
|ϕ> is the projection orbital. The projectability ∈ [0, 1] is given by |<ψ|ϕ>|^2.
The ϕ should be properly normalized, otherwise the projectability can
be larger than 1.

# Arguments
- `U`: a series of gauge matrices

# Returns
- `P`: a series of projectability vectors,
    ``p_{m k} = \\sum_{n} |U_{k m n}|^2``
    where `m` is the index of band, `k` is the index of kpoint,
    and `n` is the index of Wannier function.
"""
function compute_projectability(U::AbstractVector)
    map(U) do u
        p = u * u'
        return real(diag(p))
    end
end

"""
    $(SIGNATURES)

Sum the projectability for the given indices of orbitals.

# Arguments
- `U`: a series of gauge matrices
- `indices`: a series of indices of orbitals to be summed over, e.g.,
    `[[1, 2], [3, 4, 5, 6]]` for summing the first two orbitals into a new one,
    and the last four orbitals into another new one.

# Returns
- `P`: a series of projectability vectors, accessed by `P[ik][m, i]`
    ``p_{m i k} = \\sum_{n \\in indices[i]} |U_{k m n}|^2``
    where `m` is the index of band, `i` is the index of new orbital,
    `k` is the index of kpoint, and `n` is the index of Wannier function.
"""
function compute_projectability(U::AbstractVector, indices::AbstractVector{<:AbstractVector{<:Integer}})
    nprojs = size(U[1], 2)
    1 <= minimum(minimum.(indices)) <= maximum(maximum.(indices)) <= nprojs || error("Indices out of bounds")
    nkpts = length(U)
    nbands = size(U[1], 1)
    T = real(eltype(U[1]))
    nprojs_new = length(indices)
    P = [zeros(T, nbands, nprojs_new) for _ in 1:nkpts]
    for (ik, uk) in enumerate(U)
        for (i, idx) in enumerate(indices)
            p = uk[:, idx] * uk[:, idx]'
            P[ik][:, i] = real(diag(p))
        end
    end
    return P
end

"""Compare two structs recursively using `isapprox`."""
function isapprox_struct(a, b; kwargs...)
    for f in propertynames(a)
        va = getproperty(a, f)
        vb = getproperty(b, f)

        if va isa Vector{String}
            all(va == vb) || return false
        elseif va isa Vector
            all(isapprox.(va, vb; kwargs...)) || return false
        elseif va isa String
            va == vb || return false
        else
            isapprox(va, vb; kwargs...) || return false
        end
    end
    return true
end

"""
    $(SIGNATURES)

Check if two vectors are equivalent within a tolerance, optionally considering periodicity.

# Arguments
- `v1`: first vector
- `v2`: second vector

# Keyword Arguments
- `atol`: absolute tolerance for comparison (default: `1e-6`)
- `periodic`: whether to consider periodic boundary conditions (default: `true`)
"""
function isequiv(v1::AbstractVector, v2::AbstractVector; atol::AbstractFloat=1e-6, periodic::Bool=true)
    d = v1 - v2
    if periodic
        d -= round.(d)
    end
    return all(isapprox.(d, 0; atol))
end

"""
    $(SIGNATURES)

Create a function that compares its argument to `x` using [`isequiv`](@ref), i.e.
a function equivalent to `y -> isequiv(y, x)`.
"""
isequiv(y; kwargs...) = x -> isequiv(x, y; kwargs...)
# Cannot use the following as I want to also fix the keyword arguments.
# isequiv(x) = Base.Fix2(isequiv, x)
