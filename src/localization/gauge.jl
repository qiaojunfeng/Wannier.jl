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
- `eigenvalues`: the eigenvalues, `n_bands × n_kpoints` matrix
- `gauges`: the gauges, `n_bands × n_wannier × n_kpoints` array
"""
function transform_gauge(eigenvalues::AbstractMatrix, gauges::AbstractArray{T, 3}) where {T}
    nkpts = size(eigenvalues, 2)
    nwann = size(gauges, 2)
    H = zeros(T, nwann, nwann, nkpts)
    for ik in 1:nkpts
        U = view(gauges, :, :, ik)
        ε = view(eigenvalues, :, ik)
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
        H[:, :, ik] .= Hermitian(U' * Diagonal(ε) * U)
    end
    return H
end

"""
    $(SIGNATURES)

Transform the gauge of the (assuming a gauge-covariant) operator `O`.

For each kpoint ``\\mathbf{k}``,
```math
U^{\\dagger}_{\\mathbf{k}} O_{\\mathbf{k}} U_{\\mathbf{k}}.
```
"""
function transform_gauge(O::AbstractArray{T1, 3}, U::AbstractArray{T2, 3}) where {T1, T2}
    T = promote_type(T1, T2)
    nkpts = size(U, 3)
    nwann = size(U, 2)
    O2 = zeros(T, nwann, nwann, nkpts)
    for ik in 1:nkpts
        u = view(U, :, :, ik)
        o = view(O, :, :, ik)
        O2[:, :, ik] .= u' * o * u
    end
    return O2
end

"""
    $(SIGNATURES)

Rotate overlap `M` matrices according to gauge `U`.

For each kpoint ``\\mathbf{k}``,
```math
U_{\\mathbf{k}}^{\\dagger} M_{\\mathbf{k},\\mathbf{b}} U_{\\mathbf{k + b}}
```
"""
function transform_gauge(
        M::AbstractArray{<:Complex, 4}, kpb_k::AbstractMatrix, U::AbstractArray{<:Complex, 3}
    )
    nkpts = size(M, 4)
    nbvecs = size(M, 3)
    nwann = size(U, 2)
    T = promote_type(eltype(M), eltype(U))
    M2 = zeros(T, nwann, nwann, nbvecs, nkpts)
    for ik in 1:nkpts
        u = view(U, :, :, ik)
        for ib in 1:nbvecs
            ikpb = kpb_k[ib, ik]
            upb = view(U, :, :, ikpb)
            M2[:, :, ib, ik] .= u' * view(M, :, :, ib, ik) * upb
        end
    end
    return M2
end

"""
    $(SIGNATURES)

Rotate the gauge of a `Model`.

# Arguments
- `model`: a [`Model`](@ref)
- `U`: unitary gauge rotation matrix, `n_bands × n_wannier × n_kpoints` array

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
function rotate_gauge!(model::Model, U::AbstractArray{<:Complex, 3}; ensure_bloch_gauge::Bool = true) end

function transform_gauge(
        model::Model, U::AbstractArray{<:Complex, 3}; ensure_bloch_gauge::Bool = false
    )
    nbands = n_bands(model)
    nkpts = n_kpoints(model)
    (size(U, 1), size(U, 3)) == (nbands, nkpts) ||
        error("inconsistent size between U and model")
    nwann = size(U, 2)

    U2 = identity_gauge(eltype(U), nkpts, nwann)

    E = model.eigenvalues
    E2 = zeros(real(eltype(U)), nwann, nkpts)
    H = zeros(eltype(U), nwann, nwann)
    atol = 1.0e-8
    diag_kpts = Int[]
    for ik in 1:nkpts
        Uₖ = view(U, :, :, ik)
        H .= Uₖ' * diagm(0 => view(E, :, ik)) * Uₖ
        ϵ = diag(H)
        if norm(H - diagm(0 => ϵ)) > atol
            if ensure_bloch_gauge
                ϵ, v = eigen(H)
                U2[:, :, ik] .= v
                push!(diag_kpts, ik)
            else
                error("H is not diagonal after gauge rotation")
            end
        end
        if any(imag(ϵ) .> atol)
            error("H has non-zero imaginary part")
        end
        E2[:, ik] .= real(ϵ)
    end

    M = model.overlaps
    kpbk = model.kstencil.kpb_k
    M2 = transform_gauge(M, kpbk, U)
    if ensure_bloch_gauge && length(diag_kpts) > 0
        M2 = transform_gauge(M2, kpbk, U2)
        for ik in diag_kpts
            U2[:, :, ik] .= inv(view(U2, :, :, ik))
        end
    end

    model2 = Model(
        model.lattice,
        model.atom_positions,
        model.atom_labels,
        model.kstencil,
        M2,
        U2,
        E2,
        falses(nwann, nkpts),
    )
    return model2
end
