export transform_gauge

"""
    $(SIGNATURES)

Transform the gauge of digonal matrices (e.g., eigenvalues of Hamiltonian).

e.g., construct k-space Hamiltonian ``H(\\mathbf{k})``.
```math
H_{\\mathbf{k}} = U_{\\mathbf{k}}^\\dagger [\\epsilon_{n \\mathbf{k}}] U_{\\mathbf{k}},
```
where ``[\\epsilon_{n \\bm{k}}]`` is a diagonal matrix with
``\\epsilon_{n \\bm{k}}`` as the diagonal elements.

# Arguments
- `eigenvalues`: the eigenvalues, length-`n_kpoints` vector, each element is a
    length-`n_bands` vector
- `gauges`: the gauges, length-`n_kpoints` vector, each element is a
    `n_bands * n_wannier` matrix
"""
function transform_gauge(
    eigenvalues::AbstractVector{V}, gauges::AbstractVector{T}
) where {V<:AbstractVector,T<:AbstractMatrix}
    nkpts = length(eigenvalues)
    @assert nkpts > 0 "empty eigenvalues"
    @assert length(gauges) == nkpts "different length of eigenvalues and gauges"
    nbands = size(gauges[1], 1)
    @assert length(eigenvalues[1]) == nbands "eigenvalues have wrong n_bands"

    return map(zip(eigenvalues, gauges)) do (ε, U)
        # I need to force Hermiticity here, otherwise in some cases,
        # especially degenerate eigenvalues, the eigenvectors of Hᵏ,
        #   F = eigen(Hᵏ)
        # does not satisfy unitarity,
        #   F.vectors' ≈ F.vectors
        # and this leads to
        #   norm(F.vectors * diagm(F.values) * F.vectors' - Hᵏ) ≈ 1e-1
        # If I compute explicitly its inverse,
        #   norm(F.vectors * diagm(F.values) * inv(F.vectors) - Hᵏ) ≈ 1e-14
        # However, replacing all the `'` by `inv` is not a good idea,
        # since gauge rotation is used a lot throughout the code;
        # so I enforce Hermiticity here.
        # See also
        # https://discourse.julialang.org/t/a-b-a-is-not-hermitian-even-when-b-is/70611
        Hermitian(U' * Diagonal(ε) * U)
    end
end

"""
    $(SIGNATURES)

Transform the gauge of the (assuming a gauge-covariant) operator `O`.

For each kpoint ``\\mathbf{k}``,
```math
U^{\\dagger}_{\\mathbf{k}} O_{\\mathbf{k}} U_{\\mathbf{k}}.
```
"""
function transform_gauge(
    O::AbstractVector{V}, U::AbstractVector{T}
) where {V<:AbstractMatrix,T<:AbstractMatrix}
    nkpts = length(U)
    @assert nkpts > 0 "U must be non-empty"
    nbands, nwann = size(U[1])
    @assert length(O) == nkpts "O has wrong n_kpoints"
    @assert size(O[1], 1) == nbands "O has wrong n_bands"

    return map(zip(O, U)) do (o, u)
        u' * o * u
    end
end

"""
    $(SIGNATURES)

Rotate overlap `M` matrices according to gauge `U`.

For each kpoint ``\\mathbf{k}``,
```math
U_{\\mathbf{k}}^{\\dagger} M_{\\mathbf{k},\\mathbf{b}} U_{\\mathbf{k + b}}
```
"""
@views function transform_gauge(M::AbstractVector, kpb_k::AbstractVector, U::AbstractVector)
    nkpts = length(M)
    @assert length(M) == length(kpb_k) == length(U) "M, kpb_k, U must have same n_kpoints"
    @assert nkpts > 0 "M must be non-empty"
    nbands, nwann = size(U[1])
    nbvecs = length(M[1])
    @assert nbvecs > 0 "M must be non-empty"
    @assert (nbands, nbands) == size(M[1][1]) "M has wrong n_bands"

    return map(1:nkpts) do ik
        return map(1:nbvecs) do ib
            ikpb = kpb_k[ik][ib]
            return U[ik]' * M[ik][ib] * U[ikpb]
        end
    end
end

"""
    $(SIGNATURES)

Rotate the gauge of a `Model`.

# Arguments
- `model`: a [`Model`](@ref)
- `U`: unitary gauge rotation matrix ``U_{\\mathbf{k}}``,
    length-`n_kpoints` vector, each element is a `n_wannier * n_wannier`

# Keyword Arguments
- `ensure_bloch_gauge`:

!!! note

    This is a inplace function, meaning the dimensions of the input `model`
    should not change, thus the `n_bands` must be equal to the `n_wannier`.
    See also [`rotate_gauge`](@ref) for a non-inplace version.

    The original `model.gauges` will not be used and will be discarded;
    the `model.overlaps`, and `model.eigenvalues` will be rotated by the input `U`.
    However, since `model.eigenvalues` is not the Hamiltonian matrices but only
    their diagonal elements, the input ,
    if `diag_H = false`, this function only support rotations that keep the Hamiltonian
    in diagonal form.

    if after rotation, the Hamiltonian is not diagonal,
    then diagonalize it and save the eigenvalues to `model.eigenvalues`, and
    the inverse of the eigenvectors to `model.gauges`; otherwise, if the rotated
    Hamiltonian is not diagonal, raise error.
"""
function rotate_gauge!(model::Model, U::AbstractVector; ensure_bloch_gauge::Bool=true) end

function transform_gauge(
    model::Model, U::Vector{Matrix{T}}; ensure_bloch_gauge::Bool=false
) where {T<:Number}
    n_bands = model.n_bands
    n_kpts = model.n_kpts
    (size(U[1], 1), length(U)) == (n_bands, n_kpts) ||
        error("U must have size (n_bands, :, n_kpts)")
    # The new n_wann
    n_wann = size(U[1], 2)

    # the new gauge is just identity
    U2 = identity_U(eltype(U[1]), n_kpts, n_wann)

    # EIG
    E = model.E
    E2 = map(m -> similar(m), E)
    H = zeros(eltype(model.U[1]), n_wann, n_wann)
    # tolerance for checking Hamiltonian
    atol = 1e-8
    # all the diagonalized kpoints, used if diag_H = true
    diag_kpts = Int[]
    for ik in 1:n_kpts
        Uₖ = U[ik]
        H .= Uₖ' * diagm(0 => E[ik]) * Uₖ
        ϵ = diag(H)
        if norm(H - diagm(0 => ϵ)) > atol
            if ensure_bloch_gauge
                # diagonalize the Hamiltonian
                ϵ, v = eigen(H)
                U2[ik] = v
                push!(diag_kpts, ik)
            else
                error("H is not diagonal after gauge rotation")
            end
        end
        if any(imag(ϵ) .> atol)
            error("H has non-zero imaginary part")
        end
        E2[ik] = real(ϵ)
    end

    # MMN
    M = model.M
    kpb_k = model.bvectors.kpb_k
    M2 = transform_gauge(M, kpb_k, U)
    if ensure_bloch_gauge && length(diag_kpts) > 0
        M2 = transform_gauge(M2, kpb_k, U2)
        # the gauge matrix needs to save the inverse of the eigenvectors
        for ik in diag_kpts
            U2[ik] = inv(U2[ik])
        end
    end

    model2 = Model(
        model.lattice,
        model.atom_positions,
        model.atom_labels,
        model.kgrid,
        model.kpoints,
        model.bvectors,
        [falses(n_wann) for i in 1:n_kpts],
        M2,
        U2,
        E2,
    )
    return model2
end
