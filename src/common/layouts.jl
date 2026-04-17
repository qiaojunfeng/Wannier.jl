using Optim: Optim

export UGauge, XYGauge, WLayout, ProductLayout

"""
Layout is an abstract interface for how the optimization parameters `x` are
packed out of / into the canonical gauge array `U` of a [`Model`](@ref).

Four concrete types are used by the rewrite:

- [`UGauge`](@ref) — `x ≡ U` directly, `(n_bands, n_wannier, n_kpoints)`.
- [`XYGauge`](@ref) — disentangled form as a contiguous
  `(n_wannier² + n_bands·n_wannier) × n_kpoints` matrix.
- [`ProductLayout`](@ref) — bundle two layouts for [`SpinModel`](@ref).
- [`WLayout`](@ref) — single rotation matrix `W` in `opt_rotate`.

Every concrete `Layout` should implement the small core interface:

    encode!(x, layout, U, frozen_bands)          -> x
    decode!(U, layout, x)                        -> U
    pack_gradient!(g, layout, GU, frozen_bands)  -> g

and — once the manifold machinery lands — `manifold(layout, model, solver)`.
"""
abstract type Layout end

"""Identity layout: `x` is the gauge array `U` itself."""
struct UGauge <: Layout end

"""Disentangled layout: `x` is the packed `XY` matrix."""
struct XYGauge <: Layout end

"""Product of two layouts; used to encode `SpinModel` gauges."""
struct ProductLayout{L1 <: Layout, L2 <: Layout} <: Layout
    first::L1
    second::L2
end

"""Single `W` matrix layout for `opt_rotate`. `x ≡ W`."""
struct WLayout <: Layout end

"""
    encode!(x, layout, U, frozen_bands)

Pack the gauge array `U` into the parameter container `x` dictated by `layout`.
"""
function encode! end

"""
    decode!(U, layout, x)

Unpack the parameter container `x` back into the canonical gauge array `U`.
"""
function decode! end

"""
    pack_gradient!(g, layout, GU, frozen_bands)

Translate the canonical dΩ/dU* gradient `GU` into the layout-native gradient
buffer `g`, applying any frozen-band masking that belongs to the layout.
"""
function pack_gradient! end

# ---- UGauge --------------------------------------------------------------

encode!(x::AbstractArray{<:Complex, 3}, ::UGauge, U::AbstractArray{<:Complex, 3}, _frozen) = (copyto!(x, U); x)
decode!(U::AbstractArray{<:Complex, 3}, ::UGauge, x::AbstractArray{<:Complex, 3}) = (copyto!(U, x); U)

function pack_gradient!(
        g::AbstractArray{<:Complex, 3},
        ::UGauge,
        GU::AbstractArray{<:Complex, 3},
        frozen::AbstractMatrix{Bool},
    )
    copyto!(g, GU)
    @inbounds for ik in axes(g, 3)
        idx_f = view(frozen, :, ik)
        view(g, idx_f, :, ik) .= 0
    end
    return g
end

# ---- XYGauge -------------------------------------------------------------

function encode!(
        x::AbstractMatrix,
        ::XYGauge,
        U::AbstractArray{<:Complex, 3},
        frozen::AbstractMatrix{Bool},
    )
    X, Y = U_to_X_Y(U, frozen)
    return X_Y_to_XY!(x, X, Y)
end

function decode!(U::AbstractArray{<:Complex, 3}, ::XYGauge, x::AbstractMatrix)
    T = eltype(U)
    n_bands, n_wann, n_kpts = size(U)
    X = Array{T, 3}(undef, n_wann, n_wann, n_kpts)
    Y = Array{T, 3}(undef, n_bands, n_wann, n_kpts)
    XY_to_X_Y!(X, Y, x)
    return X_Y_to_U!(U, X, Y)
end

function pack_gradient!(
        g::AbstractMatrix,
        ::XYGauge,
        GU::AbstractArray{<:Complex, 3},
        frozen::AbstractMatrix{Bool},
    )
    # Callers that already have decoded X/Y on hand can call
    # `pack_gradient_xy!(g, GU, X, Y, frozen)` directly to avoid re-decoding.
    error(
        "pack_gradient!(g, ::XYGauge, GU, frozen) needs decoded X/Y; call " *
            "pack_gradient_xy!(g, GU, X, Y, frozen) with the stored decoded X/Y."
    )
end

"""
    pack_gradient_xy!(g, GU, X, Y, frozen)

Transform `dΩ/dU*` into the packed `XY` gradient, using the already-decoded
`X`, `Y` components. Frozen-band masking is applied inside.
"""
function pack_gradient_xy!(
        g::AbstractMatrix,
        GU::AbstractArray{<:Complex, 3},
        X::AbstractArray{<:Complex, 3},
        Y::AbstractArray{<:Complex, 3},
        frozen::AbstractMatrix{Bool},
    )
    GX, GY = GU_to_GX_GY(GU, X, Y, frozen)
    return X_Y_to_XY!(g, GX, GY)
end

# ---- ProductLayout -------------------------------------------------------

# `x` for a `ProductLayout` is a `NamedTuple{(:up, :dn), Tuple{X1, X2}}` for
# now — the concrete packing into a single contiguous buffer is handled in
# the SpinModel problem constructor (commit Q/R), not here.

encode!(x, layout::ProductLayout, (Uup, Udn), (frozen_up, frozen_dn)) = (
    encode!(x.up, layout.first,  Uup, frozen_up);
    encode!(x.dn, layout.second, Udn, frozen_dn);
    x
)

decode!((Uup, Udn), layout::ProductLayout, x) = (
    decode!(Uup, layout.first,  x.up);
    decode!(Udn, layout.second, x.dn);
    (Uup, Udn)
)

# ---- WLayout -------------------------------------------------------------

encode!(x::AbstractMatrix, ::WLayout, W::AbstractMatrix, _frozen) = (copyto!(x, W); x)
decode!(W::AbstractMatrix, ::WLayout, x::AbstractMatrix) = (copyto!(W, x); W)
pack_gradient!(g::AbstractMatrix, ::WLayout, GW::AbstractMatrix, _frozen) = (copyto!(g, GW); g)

# ---- Manifold construction ----------------------------------------------

"""
    manifold(layout, model)

Build the Optim manifold appropriate for `(layout, model)`. Solver
parameterization (e.g. `ManoptLBFGS`) is added later via a third argument;
until then the single-argument form returns an Optim.jl manifold.
"""
function manifold end

function manifold(::UGauge, model)
    nw = n_wannier(model)
    return Optim.PowerManifold(Optim.Stiefel_SVD(), (nw, nw), (n_kpoints(model),))
end

function manifold(::XYGauge, model)
    nw = n_wannier(model)
    nb = n_bands(model)
    nk = n_kpoints(model)
    per_k = Optim.ProductManifold(Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw))
    return Optim.PowerManifold(per_k, (nw^2 + nb * nw,), (nk,))
end

manifold(::WLayout, _model) = Optim.Stiefel_SVD()

function manifold(::ProductLayout{XYGauge, XYGauge}, model)
    nw = n_wannier(model.up)
    nb = n_bands(model.up)
    nk = n_kpoints(model.up)
    per_k_up = Optim.ProductManifold(Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw))
    per_k_dn = Optim.ProductManifold(Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw))
    n_inner = nw^2 + nb * nw
    k_combined = Optim.ProductManifold(per_k_up, per_k_dn, (n_inner,), (n_inner,))
    return Optim.PowerManifold(k_combined, (2 * n_inner,), (nk,))
end

# ---- Legacy conversion helpers ------------------------------------------

"""
    X_Y_to_U(X, Y)

Convert the `(X, Y)` layout to the `U` layout.

There are three formats: `U`, `(X, Y)`, and `XY` stored contiguously in memory.
For each kpoint,
- `U`: `size(U) = (n_bands, n_wann, n_kpts)`, the format used in the rest of the code
- `(X, Y)`: `size(X) = (n_wann, n_wann, n_kpts)`, `size(Y) = (n_bands, n_wann, n_kpts)`, intermediate format
- `XY`: this is the format used in the optimizer
"""
function X_Y_to_U(X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T <: Complex}
    n_bands = size(Y, 1)
    n_wann = size(Y, 2)
    n_kpts = size(Y, 3)

    U = zeros(T, n_bands, n_wann, n_kpts)
    return X_Y_to_U!(U, X, Y)
end

function X_Y_to_U!(U::AbstractArray{T, 3}, X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T}
    @inbounds for ik in axes(U, 3)
        mul!(view(U, :, :, ik), view(Y, :, :, ik), view(X, :, :, ik))
    end
    return U
end

"""
    U_to_X_Y(U, frozen)

Convert the `U` layout to the `(X, Y)` layout.

See also [`X_Y_to_U`](@ref).

# Arguments
- `U`: `n_bands × n_wann × n_kpts`
- `frozen`: `n_bands × n_kpts` BitMatrix
"""
function U_to_X_Y(U::AbstractArray{T, 3}, frozen::AbstractMatrix{Bool}) where {T <: Complex}
    nkpts = size(U, 3)
    nbands, nwann = size(U, 1), size(U, 2)

    X = zeros(T, nwann, nwann, nkpts)
    Y = zeros(T, nbands, nwann, nkpts)

    @inbounds for ik in 1:nkpts
        idx_f = view(frozen, :, ik)
        idx_nf = .!idx_f
        n_froz = count(idx_f)

        Af = orthonorm_freeze(view(U, :, :, ik), idx_f)
        Uf = Af[idx_f, :]
        Ur = Af[idx_nf, :]

        # determine Y
        Y[idx_f, 1:n_froz, ik] .= Matrix{T}(I, n_froz, n_froz)

        if n_froz != nwann
            Pr = Ur * Ur'
            Pr = Hermitian((Pr + Pr') / 2)
            D, V = eigen(Pr) # sorted by increasing eigenvalue
            Y[idx_nf, (n_froz + 1):end, ik] .= V[:, (end - nwann + n_froz + 1):end]
        end

        # determine X
        X[:, :, ik] .= orthonorm_lowdin(view(Y, :, :, ik)' * Af)
    end

    return X, Y
end

"""
    XY_to_X_Y(XY, n_bands, n_wann)

Convert the `XY` layout to the `(X, Y)` layout.

See also [`X_Y_to_U`](@ref).

# Arguments
- `XY`: `n_bands * n_wann * n_kpts` contiguous array
- `n_bands`: number of bands, to be used to reshape `XY`
- `n_wann`: number of wannier functions, to be used to reshape `XY`
"""
function XY_to_X_Y(XY::AbstractMatrix{T}, nbands::Int, nwann::Int) where {T <: Complex}
    nkpts = size(XY, 2)

    X = zeros(T, nwann, nwann, nkpts)
    Y = zeros(T, nbands, nwann, nkpts)
    return XY_to_X_Y!(X, Y, XY)
end

function XY_to_X_Y!(X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}, XY::AbstractMatrix) where {T}
    n_wann2 = size(X, 1)^2
    @inbounds for ik in axes(X, 3)
        xk = view(X, :, :, ik)
        yk = view(Y, :, :, ik)
        for i in eachindex(xk)
            xk[i] = XY[i, ik]
        end
        for i in eachindex(yk)
            yk[i] = XY[n_wann2 + i, ik]
        end
    end
    return X, Y
end

"""
    X_Y_to_XY(X, Y)

Convert the `(X, Y)` layout to the `XY` layout.

See also [`X_Y_to_U`](@ref).
"""
function X_Y_to_XY(X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T <: Complex}
    nkpts = size(Y, 3)
    nbands, nwann = size(Y, 1), size(Y, 2)
    n = nwann^2
    XY = zeros(T, n + nbands * nwann, nkpts)
    return X_Y_to_XY!(XY, X, Y)
end

function X_Y_to_XY!(XY::AbstractMatrix, X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T}
    @inbounds for ik in axes(X, 3)
        xk = view(X, :, :, ik)
        yk = view(Y, :, :, ik)
        n = length(xk)
        for i in eachindex(xk)
            XY[i, ik] = xk[i]
        end
        for i in eachindex(yk)
            XY[n + i, ik] = yk[i]
        end
    end
    return XY
end

@doc raw"""
    GU_to_GX_GY(G, X, Y, frozen)

Compute dΩ/dX and dΩ/dY from dΩ/dU.

Acutally they are the conjugate gradients, e.g., ``\frac{d \Omega}{d U^*}``.

# Arguments
- `G`: `n_bands × n_wann × n_kpts` array for gradient dΩ/dU
- `X`: `n_wann × n_wann × n_kpts` array for X
- `Y`: `n_bands × n_wann × n_kpts` array for Y
- `frozen`: `n_bands × n_kpts` BitMatrix for frozen bands
"""
function GU_to_GX_GY(
        G::AbstractArray{T, 3},
        X::AbstractArray{T, 3},
        Y::AbstractArray{T, 3},
        frozen::AbstractMatrix{Bool},
    ) where {T}
    n_kpts = size(X, 3)
    GX = zeros(T, size(X)...)
    GY = zeros(T, size(Y)...)

    @inbounds for ik in 1:n_kpts
        idx_f = view(frozen, :, ik)
        n_froz = count(idx_f)

        mul!(view(GX, :, :, ik), view(Y, :, :, ik)', view(G, :, :, ik))
        mul!(view(GY, :, :, ik), view(G, :, :, ik), view(X, :, :, ik)')

        view(GY, idx_f, :, ik) .= 0
        view(GY, :, 1:n_froz, ik) .= 0
    end

    return GX, GY
end
