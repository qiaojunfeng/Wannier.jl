using LinearAlgebra
using Optim: Optim

"""
    get_frozen_bands(E, dis_froz_max, dis_froz_min)

Generate a `BitMatrix` of frozen bands by checking the two frozen windows.

# Arguments
- `E`: the energy eigenvalues, `n_bands × n_kpoints` matrix
- `dis_froz_max`: the upper bound of the frozen window
- `dis_froz_min`: the lower bound of the frozen window

!!! note

    The `dis_froz_max` and `dis_froz_min` work similarly as `Wannier90`.
"""
function get_frozen_bands(E::AbstractMatrix, dis_froz_max, dis_froz_min = -Inf)
    return BitMatrix((E .>= dis_froz_min) .& (E .<= dis_froz_max))
end

#TODO I don't think this works; degen not defined
"""
    set_frozen_degen!(frozen_bands, E, atol=1e-4)

Freeze bands which are degenerate with the `frozen_bands`.

In some cases, we might want to freeze the whole set of degenerated eigen vectors.

# Arguments
- `frozen_bands`: the `BitMatrix` of frozen bands, `n_bands × n_kpoints`
- `E`: the energy eigenvalues, `n_bands × n_kpoints` matrix
- `atol`: the tolerance of degeneracy
"""
function set_frozen_degen!(
        frozen_bands::AbstractMatrix{Bool}, E::AbstractMatrix, atol::T = 1.0e-4
    ) where {T <: Real}
    nbands, n_kpts = size(E)
    atol <= 0 && error("atol must be positive")

    for ik in 1:n_kpts
        frozen_k = view(frozen_bands, :, ik)

        # if cluster of eigenvalues and count(frozen_k) > 0, take them all
        if degen && count(frozen_k) > 0
            ib = findlast(frozen_k)

            while ib < nbands
                if E[ib + 1, ik] < E[ib, ik] + atol
                    ib += 1
                    frozen_k[ib] = true
                else
                    break
                end
            end
        end
    end

    return nothing
end

"""
    check_frozen_bands(frozen_bands, n_wann)

Sanity check that the number of frozen bands at each kpoint <= `n_wann`.

# Arguments
- `frozen_bands`: the `BitMatrix` of frozen bands, `n_bands × n_kpoints`
- `n_wann`: the number of wannier functions
"""
function check_frozen_bands(frozen_bands::AbstractMatrix{Bool}, n_wann::Int)
    nbands, n_kpts = size(frozen_bands)
    n_wann > nbands && error("n_wann > nbands")

    for ik in 1:n_kpts
        if count(view(frozen_bands, :, ik)) > n_wann
            error("Too many frozen bands")
        end
    end
    return
end

"""
    set_frozen_win!(model, dis_froz_max, dis_froz_min=-Inf; degen=false, degen_atol=1e-4)

Set frozen bands of the `Model` according to two energy windows.

# Arguments
- `model`: the `Model` to be set
- `dis_froz_max`: the upper bound of the frozen window
- `dis_froz_min`: the lower bound of the frozen window

# Keyword Arguments
- `degen`: whether to freeze the whole set of degenerated eigen vectors
- `degen_atol`: the tolerance of degeneracy
"""
function set_frozen_win!(
        model::Model{T},
        dis_froz_max::T,
        dis_froz_min::T = -Inf;
        degen::Bool = false,
        degen_atol::T = 1.0e-4,
    ) where {T <: Real}
    frozen_bands = get_frozen_bands(model.eigenvalues, dis_froz_max, dis_froz_min)

    if degen
        set_frozen_degen!(frozen_bands, model.eigenvalues, degen_atol)
    end

    check_frozen_bands(frozen_bands, n_wannier(model))

    model.frozen_bands .= frozen_bands

    return nothing
end

"""
    get_frozen_proj(E, U, dis_proj_max)

Get frozen bands according to band projectability.

# Arguments
- `E`: the energy eigenvalues, `n_bands × n_kpoints` matrix
- `U`: the gauge rotation matrices, `n_bands × n_wann × n_kpoints` array
- `dis_proj_max`: the upper bound projectability.
    Bands with projectability >= `dis_proj_max` are frozen.

!!! note

    The band projectability for band ``n`` at kpoint ``\\bm{k}`` is calculated by
    ``p_{n \\bm{k}} = \\sum_{m=1}^{m=n_{wann}} | U_{nm \\bm{k}} |^2``,
    and usually each element ``p_{n \\bm{k}} \\in [0.0, 1.0]``.
    In such cases, the `dis_proj_max` is usually set to sth. like `0.9` to
    freeze high-projectability bands.
"""
function get_frozen_proj(
        E::AbstractMatrix{T},
        U::AbstractArray{Complex{T}, 3},
        dis_proj_max::T,
    ) where {T <: Real}
    nbands, n_kpts = size(E)
    frozen_bands = falses(nbands, n_kpts)

    for ik in 1:n_kpts
        Uₖ = view(U, :, :, ik)
        p = dropdims(real(sum(Uₖ .* conj(Uₖ); dims = 2)); dims = 2)
        frozen_bands[p .>= dis_proj_max, ik] .= true
    end

    return frozen_bands
end

"""
    set_frozen_proj!(model, dis_proj_max; degen=false, degen_atol=1e-4)

Set frozen bands of the `Model` according to projectability.

# Arguments
- `model`: the `Model` to be set
- `dis_proj_max`: the upper bound projectability.
    Bands with projectability >= `dis_proj_max` are frozen.

# Keyword Arguments
- `degen`: whether to freeze the whole set of degenerated eigen vectors
- `degen_atol`: the tolerance of degeneracy
"""
function set_frozen_proj!(
        model::Model{T}, dis_proj_max::T; degen::Bool = false, degen_atol::T = 1.0e-4
    ) where {T <: Real}
    frozen_bands = get_frozen_proj(model.eigenvalues, model.gauges, dis_proj_max)

    if degen
        set_frozen_degen!(frozen_bands, model.eigenvalues, degen_atol)
    end

    check_frozen_bands(frozen_bands, n_wannier(model))

    model.frozen_bands .= frozen_bands

    return nothing
end

"""
    orthonormalize_frozen(U, frozen)

Normalize and freeze a block of a matrix.

Conditions:
- Block form:  `U = vcat(Uf, Ur)`
- Semiunitary
  - `U' * U = I`
  - `Uf' * Uf + Ur' * Ur = I`
- Frozen:      `Uf * Uf' = I`
- Also:        `Uf * Ur' = 0`

Strategy:
1. orthogonalize `Uf`
2. project `Uf` out of `Ur`
3. orthogonalize the range of `Ur`

# Arguments
- `U`: the matrix to be orthonormalized and frozen
- `frozen`: the `BitVector` specifying which bands are frozen
"""
function orthonormalize_frozen(U::AbstractMatrix{T}, frozen::AbstractVector{Bool}) where {T <: Complex}
    nbands, nwann = size(U)
    non_frozen = .!frozen

    # Make sure Uf can fully represent frozen bands.
    # Uf = <ψ|g>, where |ψ> is Bloch wfcs, |g> is the guiding functions.
    # We do a Lowdin orthonormalization on Uf so that Uf * Uf' = I,
    # i.e. <ψ|g'><g'|ψ> = I -> |g'>s span the frozen |ψ>s.
    Uf = U[frozen, :]
    Uf = lowdin_orthonormalize(Uf)

    # Remove Uf out of Ur, i.e. do not destroy frozen space
    # The projector of the frozen states represented on the |g> basis is
    #     |ψf><ψf| = |g><g|ψf><ψf|g><g| = |g> Uf' * Uf <g|
    # The projector of the non-frozen states on GF basis is
    #     |ψr><ψr| = |g><g|ψr><ψr|g><g| = |g> Ur' * Ur <g|
    # The space |ψr><ψr| should not destroy the frozen space |ψf><ψf|.
    # To achieve this, we remove |ψf><ψf| components out of |ψr><ψr|
    #     |ψr> -= |ψf><ψf|ψr>
    # =>  |g><g|ψr> -= |g><g|ψf><ψf|g><g|ψr>
    # =>  Ur' -= Uf' * Uf * Ur'
    Ur = U[non_frozen, :]
    Ur -= Ur * Uf' * Uf

    # alternative method, maybe more stable but slower
    # M = Ur' * Ur
    # S, U = eig(Hermitian(M))
    # S = real(S)
    # @assert !any(x -> 1e-11 <= x <= 1e-9, S)
    # @assert count(x -> x > 1e-10, S) == size(A,2) - nfrozen
    # Sm12 = map(x-> x < 1e-10 ? 0. : 1/sqrt(x), S)
    # Ur = Ur*(U*diagm(Sm12)*U')

    # Renormalize the range of Ur
    # The remaining Wannier function (WF) |wr> = |ψ> Ur,
    # after removal of Uf, we need to renormalize so the |wr> are orthonormal.
    # I = <wr|wr> = Ur' <ψ|ψ> Ur  =>  Ur' * Ur = I
    # Use Lowdin normalization but needs to limit the number of independent vectors.
    A, S, B = svd(Ur)
    atol = 1.0e-10
    @assert count(x -> x > atol, S) == nwann - count(frozen) "S = $S"
    S[S .> atol] .= 1
    S[S .< atol] .= 0
    Ur = A * Diagonal(S) * B'

    V = similar(U)
    V[frozen, :] .= Uf
    V[non_frozen, :] .= Ur

    # Semiunitary
    @assert isapprox(V' * V, I; atol)
    # Frozen
    @assert isapprox(V[frozen, :] * V[frozen, :]', I; atol)
    # Independent
    @assert norm(Uf * Ur') < atol

    return V
end

"""
    zero_froz_grad!(G, frozen)

Set gradient of frozen bands to 0.

This is used in test.

# Arguments
- `G`: gradient of the spread, in `XY` layout
- `frozen`: `BitMatrix` for frozen bands, `n_bands × n_kpts`
"""
function zero_froz_grad!(G::AbstractMatrix, frozen::AbstractMatrix{Bool})
    nbands, nkpts = size(frozen)
    size(G, 2) == nkpts || error("length(G) != n_kpts")
    nwann = round(Int, (-nbands + sqrt(nbands^2 + 4 * size(G, 1))) / 2)

    GX, GY = Wannier.XY_to_X_Y(G, nbands, nwann)
    @inbounds for ik in 1:nkpts
        idx_f = view(frozen, :, ik)
        n_froz = count(idx_f)
        view(GY, idx_f, :, ik) .= 0
        view(GY, :, 1:n_froz, ik) .= 0
    end
    G .= Wannier.X_Y_to_XY(GX, GY)
    return nothing
end

