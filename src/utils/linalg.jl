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

function orthonorm_lowdin(U::AbstractArray{T, 3}) where {T}
    U2 = similar(U)
    for ik in axes(U, 3)
        U2[:, :, ik] .= orthonorm_lowdin(view(U, :, :, ik))
    end
    return U2
end

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
function powm(U::AbstractMatrix{T}, p::F) where {T <: Union{Complex, Real}, F <: Real}
    # Workaround, eigen incompatible with lazy adjoint.
    d, V = eigen(Matrix(U))

    V = orthonorm_lowdin(V)
    # accuracy = norm(V * Diagonal(d) * V' - U)
    # @assert accuracy < 1e-10

    return V * Diagonal(d .^ p) * V'
end

"""
    eig_log(O; guiding=false, tol=0.01)

Eigendecomposition-based matrix logarithm of a (near-)unitary matrix `O`.

Return `(V, logd)` such that `log(O) ≈ V * Diagonal(logd) * V'`, with `V`
Löwdin-orthonormalized so the reconstruction is robust to degenerate
eigenvalues. Since the eigenvector column phases cancel in `V * D * V'`, this
matches the raw `eigen` result in the non-degenerate case.

Branch of the logarithm (which `2π` window each eigenphase lands in):

- Default (`guiding=false`): principal branch `angle(d) ∈ (-π, π]`, but
  eigenphases within `tol` of `-π` are lifted by `+2π`. This keeps a cluster of
  eigenvalues near `-1` on the same side of the cut and reproduces the legacy
  heuristic used by [`parallel_transport`](@ref).
- `guiding=true`: place the cut in the *largest empty arc* of the eigenvalue
  distribution (the intrinsic analog of [`choose_pole`](@ref) picking a pole far
  from the data). NOTE: this is *not* generally better — the branch that yields
  the smoothest gauge depends on the propagation history, not on the spectrum of
  `O` alone, so this option can degrade the result. Kept for experimentation.
"""
function eig_log(O::AbstractMatrix; guiding::Bool = false, tol::Real = 0.01)
    d, V = eigen(Matrix(O))
    V = orthonorm_lowdin(V)
    θ = angle.(d)
    if guiding && length(θ) > 1
        θs = sort(θ)
        n = length(θs)
        # start with the wrap-around gap, then scan the interior gaps
        gmax = θs[1] + 2π - θs[n]
        cut = θs[n] + gmax / 2
        for i in 1:(n - 1)
            g = θs[i + 1] - θs[i]
            if g > gmax
                gmax = g
                cut = θs[i] + g / 2
            end
        end
        # lift each eigenphase into (cut - 2π, cut], i.e. away from the cut
        θ = @. cut - mod(cut - θ, 2π)
    else
        # principal branch, keeping near-(-π) eigenphases on the +π side
        θ = @. θ + 2π * (θ < -π + tol)
    end
    return V, im .* θ
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
    U = zeros(T, nwann, nwann, nkpts)
    for ik in 1:nkpts
        for i in 1:nwann
            U[i, i, ik] = one(T)
        end
    end
    return U
end

"""
    zeros_gauge(T, nkpts, nbands, nwann)
    zeros_gauge(T, nkpts, nwann)

Allocate gauge matrix `U` filled with zeros.

See also [`identity_gauge`](@ref).
"""
function zeros_gauge end

function zeros_gauge(T::Type, nkpts::Integer, nbands::Integer, nwann::Integer)
    return zeros(T, nbands, nwann, nkpts)
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
function merge_gauge(U::AbstractArray{T1, 3}, V::AbstractArray{T2, 3}) where {T1, T2}
    T = promote_type(T1, T2)
    nkpts = size(U, 3)
    nbands = size(U, 1)
    nwann = size(V, 2)
    W = zeros(T, nbands, nwann, nkpts)
    for ik in 1:nkpts
        W[:, :, ik] .= view(U, :, :, ik) * view(V, :, :, ik)
    end
    return W
end

"""Apply the same `W` rotation to every k-point of `U`."""
function merge_gauge(U::AbstractArray{T1, 3}, W::AbstractMatrix{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    nbands, nwann_in, nkpts = size(U)
    size(W, 1) == nwann_in || error("W rows ($(size(W, 1))) must equal n_wann_in ($nwann_in)")
    nwann_out = size(W, 2)
    UW = zeros(T, nbands, nwann_out, nkpts)
    for ik in 1:nkpts
        UW[:, :, ik] .= view(U, :, :, ik) * W
    end
    return UW
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
    U = zeros(T, nbands, nwann, nkpts)
    for ik in 1:nkpts
        U[:, :, ik] .= rand_gauge(T, nbands, nwann)
    end
    return U
end

"""
    $(SIGNATURES)

Allocate overlap `M` matrices filled with zeros.

The `M` can be accessed by `view(M,:,:,ib,ik)`, where `ik`, `ib`, `m`, `n` are
the indices of kpoints, b-vectors, bands, and WFs, respectively.
"""
function zeros_overlap(T::Type, nkpts::Integer, nbvecs::Integer, nbands::Integer)
    return zeros(T, nbands, nbands, nbvecs, nkpts)
end

"""
    $(SIGNATURES)

Allocate eigenvalues matrix filled with zeros.

The returned `E` can be accessed by `E[n, ik]`, where `n`, `ik` are the indices
of bands and kpoints, respectively.
"""
function zeros_eigenvalues(T::Type, nkpts::Integer, nwann::Integer)
    return zeros(T, nwann, nkpts)
end

"""
    $(SIGNATURES)

Check if matrix is unitary or semi-unitary for all the kpoints.

I.e. if it has orthogonal columns.
"""
function isunitary(U::AbstractArray{T, 3}; atol = 1.0e-10) where {T}
    for ik in axes(U, 3)
        u = view(U, :, :, ik)
        if norm(u' * u - I) > atol
            @debug "not unitary" ik norm(u' * u - I)
            return false
        end
    end
    return true
end

"""
    $(SIGNATURES)

The total projectability (summed over all the orbitals) at each kpoint.

The projections (gauge matrices) are <ψ|ϕ> where |ψ> is the Bloch state and
|ϕ> is the projection orbital. The projectability ∈ [0, 1] is given by |<ψ|ϕ>|^2.
The ϕ should be properly normalized, otherwise the projectability can
be larger than 1.

# Arguments
- `U`: gauge matrices, indexed by `U[k][m, n]` where `k`, `m`, `n` are the
    indices of kpoints, bands, and WFs, respectively.

# Returns
- `P`: total projectabilities ``p_{m k} = \\sum_{n} |U_{k m n}|^2``
"""
function total_projectability(U::AbstractVector)
    return map(U) do u
        p = u * u'
        return real(diag(p))
    end
end

"""
    $(SIGNATURES)

Sum the projectability for the given indices of orbitals.

# Arguments
- `U`: gauge matrices, indexed by `U[k][m, n]` where `k`, `m`, `n` are the
    indices of kpoints, bands, and WFs, respectively.
- `index_map`: a series of indices of orbitals to be summed over, e.g.,
    `[[1, 2], [3, 4, 5, 6]]` for summing the first two orbitals into a new one,
    and the last four orbitals into another new one.

# Returns
- `P`: orbital-resolved projectabilities, indexed by `P[k][m, i]`
    ``p_{m i k} = \\sum_{n \\in \\text{index\\_map}[i]} |U_{k m n}|^2``
    where `m` is the index of band, `i` is the index of new orbital,
    `k` is the index of kpoint, and `n` is the index of Wannier function.
"""
function projectability(U::AbstractVector, index_map::AbstractVector{<:AbstractVector{<:Integer}})
    nprojs = size(U[1], 2)
    1 <= minimum(minimum.(index_map)) <= maximum(maximum.(index_map)) <= nprojs || error("Indices out of bounds")
    nkpts = length(U)
    nbands = size(U[1], 1)
    T = real(eltype(U[1]))
    nprojs_new = length(index_map)
    P = [zeros(T, nbands, nprojs_new) for _ in 1:nkpts]
    for (ik, uk) in enumerate(U)
        for (i, idx) in enumerate(index_map)
            p = uk[:, idx] * uk[:, idx]'
            P[ik][:, i] = real(diag(p))
        end
    end
    return P
end

"""
Projectability for each eigenstate onto each orbital, accessed by P[ip][ib, ik].
"""
function projectability(U::AbstractVector)
    index_map = [[i] for i in 1:size(U[1], 2)]
    return projectability(U, index_map)
end

function projectability(U::AbstractArray{<:Any, 3})
    return projectability([view(U, :, :, ik) for ik in axes(U, 3)])
end

function projectability(U::AbstractArray{<:Any, 3}, index_map::AbstractVector{<:AbstractVector{<:Integer}})
    return projectability([view(U, :, :, ik) for ik in axes(U, 3)], index_map)
end

"""
    $(SIGNATURES)

Sum the projectability over unique labels.

# Arguments
- `P`: the projectability indexed by `P[k][m, i]` where `k`, `m`, `i` are the indices of
    kpoints, bands, and orbitals, respectively. Note this is NOT the gauge matrices, they
    should be real-valued matrices.
- `labels`: the labels of the orbitals, indexed by `labels[i]` where `i` is the
    index of orbital. Orbitals with the same label will be summed together.
    e.g., `["s", "p", "p", "p"]` for keeping the first orbital as "s", and summing the last three orbitals into "p".
    One can also use `["Si1:s", "Si1:p", "Si1:p", "Si1:p", "Si2:s", "Si2:p", "Si2:p", "Si2:p"]`
    for summing the orbitals of same angular momentum on the same atom, but
    keeping the orbitals of different atoms separate.
"""
function sum_projectability(P::AbstractVector, labels::AbstractVector)
    nkpts = length(P)
    nbands, nprojs = size(P[1])
    labels_new = unique(labels)
    nprojs_new = length(labels_new)
    T = eltype(P[1])
    P_new = [zeros(T, nbands, nprojs_new) for _ in 1:nkpts]
    orb_idxs = map(labels_new) do lab
        findall(==(lab), labels)
    end
    for (ik, pk) in enumerate(P)
        for i in 1:nprojs_new
            idx = orb_idxs[i]
            P_new[ik][:, i] = sum(pk[:, idx]; dims = 2)
        end
    end
    return P_new, labels_new
end

"""
Sum projectability over all orbitals.
"""
function sum_projectability(P::AbstractVector)
    nprojs = size(P[1], 2)
    labels = ["" for _ in 1:nprojs]
    P = sum_projectability(P, labels)[1]
    return [p[:, 1] for p in P]
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
function isequiv(v1::AbstractVector, v2::AbstractVector; atol::AbstractFloat = 1.0e-6, periodic::Bool = true)
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
