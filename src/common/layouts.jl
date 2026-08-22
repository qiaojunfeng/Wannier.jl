using Optim: Optim

export ULayout, XYLayout, WLayout, ProductLayout

"""
Layout is an abstract interface for how the optimization parameters `x` are
packed out of / into the canonical gauge array `U` of a [`Model`](@ref).

Four concrete types are used by the rewrite:

- [`ULayout`](@ref) — `x ≡ U` directly, `(n_bands, n_wannier, n_kpoints)`.
- [`XYLayout`](@ref) — disentangled form as a contiguous
  `(n_wannier² + n_bands·n_wannier) × n_kpoints` matrix.
- [`ProductLayout`](@ref) — bundle two layouts for [`SpinModel`](@ref).
- [`WLayout`](@ref) — single rotation matrix `W` in `opt_rotate`.

Every concrete `Layout` should implement the small core interface:

    encode!(x, layout, U, frozen_bands)          -> x
    decode!(U, layout, x)                        -> U
    encode_gradient!(g, layout, GU, frozen_bands)  -> g

and — once the manifold machinery lands — `manifold(layout, model, solver)`.
"""
abstract type Layout end

"""Identity layout: `x` is the gauge array `U` itself."""
struct ULayout <: Layout end

"""Disentangled layout: `x` is the packed `XY` matrix."""
struct XYLayout <: Layout end

"""Product of two layouts; used to encode `SpinModel` gauges."""
struct ProductLayout{L1 <: Layout, L2 <: Layout} <: Layout
    first::L1
    second::L2
end

"""Single `W` matrix layout for `opt_rotate`. `x ≡ W`."""
struct WLayout <: Layout end

"""
    initial_x(layout, model) -> x

Build the starting parameter container for `layout` out of `model.gauges`.
"""
function initial_x end

"""
    decode!(layout, x, model, ws) -> U

Convert the layout-native parameters `x` into the canonical gauge that
objectives consume: an `(n_bands, n_wannier, n_kpoints)` array, or the
`(up, dn)` pair for a [`ProductLayout`](@ref). Intermediates that the gradient
encoding needs again — the `X`/`Y` blocks, the decoded `U` — are stashed in
`ws`, so `encode_gradient!` does not recompute them.
"""
function decode! end

"""
    encode_gradient!(g, layout, model, ws) -> g

Translate the canonical gradient `dΩ/dU*`, which the objective left in `ws.GU`,
into the layout-native gradient buffer `g`. Frozen-band masking belongs here —
no other part of the code applies it.
"""
function encode_gradient! end

# ---- ULayout --------------------------------------------------------------

initial_x(::ULayout, model) = copy(model.gauges)

# x ≡ U, so decoding is the identity — no copy, the objective reads `x` directly.
decode!(::ULayout, x::AbstractArray{<:Complex, 3}, model, ws) = x

encode_gradient!(g::AbstractArray{<:Complex, 3}, ::ULayout, model, ws) = copyto!(g, ws.GU)

# ---- XYLayout -------------------------------------------------------------

function initial_x(::XYLayout, model)
    X, Y = U_to_X_Y(model.gauges, model.frozen_bands)
    return X_Y_to_XY(X, Y)
end

function decode!(::XYLayout, x::AbstractMatrix, model, ws)
    X, Y = XY_to_X_Y!(ws.X, ws.Y, x)
    return X_Y_to_U!(ws.U, X, Y)
end

encode_gradient!(g::AbstractMatrix, ::XYLayout, model, ws) =
    encode_gradient_xy!(g, ws.GU, ws.X, ws.Y, model.frozen_bands)

"""
    encode_gradient_xy!(g, GU, X, Y, frozen)

Transform `dΩ/dU*` into the packed `XY` gradient, using the already-decoded
`X`, `Y` blocks. Frozen-band masking is applied inside.
"""
function encode_gradient_xy!(
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

# `x` stacks the two channels' packed XY blocks into one `(2 n_inner, n_kpoints)`
# matrix, up first. `_inner_size` is the per-channel row count.
_inner_size(model) = n_wannier(model)^2 + n_bands(model) * n_wannier(model)

function initial_x(::ProductLayout, model)
    Xup, Yup = U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn, Ydn = U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
    return vcat(X_Y_to_XY(Xup, Yup), X_Y_to_XY(Xdn, Ydn))
end

function decode!(::ProductLayout, x, model, ws)
    ni = _inner_size(model.up)
    xr = reshape(x, (2 * ni, n_kpoints(model)))
    Xup, Yup = XY_to_X_Y!(ws.up.X, ws.up.Y, @view xr[1:ni, :])
    Xdn, Ydn = XY_to_X_Y!(ws.dn.X, ws.dn.Y, @view xr[(ni + 1):end, :])
    return (X_Y_to_U!(ws.up.U, Xup, Yup), X_Y_to_U!(ws.dn.U, Xdn, Ydn))
end

function encode_gradient!(g, ::ProductLayout, model, ws)
    ni = _inner_size(model.up)
    gr = reshape(g, (2 * ni, n_kpoints(model)))
    encode_gradient_xy!(
        view(gr, 1:ni, :), ws.up.GU, ws.up.X, ws.up.Y, model.up.frozen_bands
    )
    encode_gradient_xy!(
        view(gr, (ni + 1):(2 * ni), :), ws.dn.GU, ws.dn.X, ws.dn.Y, model.dn.frozen_bands
    )
    return g
end

# ---- WLayout -------------------------------------------------------------

# A single rotation shared by all kpoints; the starting point is the identity.
initial_x(::WLayout, model) =
    Matrix{eltype(model.gauges)}(I, n_wannier(model), n_wannier(model))

# The canonical gauge is the same rotation applied at every kpoint: UW_k = U_k W.
function decode!(::WLayout, W::AbstractMatrix, model, ws)
    @inbounds for ik in axes(ws.U, 3)
        mul!(view(ws.U, :, :, ik), view(model.gauges, :, :, ik), W)
    end
    return ws.U
end

# Chain rule for UW_k = U_k W collapses the k index: dΩ/dW* = Σ_k U_k† dΩ/dUW_k*.
function encode_gradient!(g::AbstractMatrix, ::WLayout, model, ws)
    fill!(g, 0)
    @inbounds for ik in axes(ws.GU, 3)
        mul!(g, view(model.gauges, :, :, ik)', view(ws.GU, :, :, ik), true, true)
    end
    return g
end

"""
    decode(layout, x, model) -> U

Non-mutating counterpart of [`decode!`](@ref): turn final optimizer output into
freshly allocated canonical gauges, with no aliasing of workspace buffers.
"""
function decode end

decode(::ULayout, x, model) = x
decode(::WLayout, x, model) = x

function decode(::XYLayout, x, model)
    X, Y = XY_to_X_Y(x, n_bands(model), n_wannier(model))
    return X_Y_to_U(X, Y)
end

function decode(::ProductLayout, x, model)
    ni = _inner_size(model.up)
    nb, nw = n_bands(model.up), n_wannier(model.up)
    Xup, Yup = XY_to_X_Y(x[1:ni, :], nb, nw)
    Xdn, Ydn = XY_to_X_Y(x[(ni + 1):end, :], nb, nw)
    return X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn)
end

# ---- Manifold construction ----------------------------------------------

"""
    manifold(layout, model)

Build the Optim manifold appropriate for `(layout, model)`. Solver
parameterization (e.g. `ManoptLBFGS`) is added later via a third argument;
until then the single-argument form returns an Optim.jl manifold.
"""
function manifold end

function manifold(::ULayout, model)
    nw = n_wannier(model)
    return Optim.PowerManifold(Optim.Stiefel_SVD(), (nw, nw), (n_kpoints(model),))
end

function manifold(::XYLayout, model)
    nw = n_wannier(model)
    nb = n_bands(model)
    nk = n_kpoints(model)
    per_k = Optim.ProductManifold(Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw))
    return Optim.PowerManifold(per_k, (nw^2 + nb * nw,), (nk,))
end

manifold(::WLayout, _model) = Optim.Stiefel_SVD()

function manifold(::ProductLayout{XYLayout, XYLayout}, model)
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

        Af = orthonormalize_frozen(view(U, :, :, ik), idx_f)
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
        X[:, :, ik] .= lowdin_orthonormalize(view(Y, :, :, ik)' * Af)
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
