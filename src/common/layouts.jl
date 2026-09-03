using Optim: Optim

export ULayout, XYLayout, WLayout, ProductLayout

"""
Layout is an abstract interface for parameterizing the canonical gauge array
`U` of a [`Model`](@ref) by the optimizer variables `x`.

Four concrete types are used by the rewrite:

- [`ULayout`](@ref) — `x ≡ U` directly, `(n_bands, n_wannier, n_kpoints)`.
- [`XYLayout`](@ref) — disentangled form as one contiguous vector containing
  each full `X` and only the active, nonfrozen part of each `Y`.
- [`ProductLayout`](@ref) — bundle two layouts for [`SpinModel`](@ref).
- [`WLayout`](@ref) — single rotation matrix `W` in `opt_rotate`.

Every concrete `Layout` implements the small core interface
[`initial_parameters`](@ref), [`assemble_gauge!`](@ref),
[`pullback_gradient!`](@ref), [`finalize_result`](@ref), and
[`manifold`](@ref).
"""
abstract type Layout end

"""Identity layout: `x` is the gauge array `U` itself."""
struct ULayout <: Layout end

"""Disentangled layout: `x` compactly packs the active `X`/`Y` blocks."""
struct XYLayout <: Layout end

"""Product of two layouts; used to parameterize `SpinModel` gauges."""
struct ProductLayout{L1 <: Layout, L2 <: Layout} <: Layout
    first::L1
    second::L2
end

"""Single `W` matrix layout for `opt_rotate`. `x ≡ W`."""
struct WLayout <: Layout end

"""
    initial_parameters(layout, model) -> x

Build the starting parameter container for `layout` out of `model.gauges`.
"""
function initial_parameters end

"""
    assemble_gauge!(layout, x, model, ws) -> U

Convert the layout-native parameters `x` into the canonical gauge that
objectives consume: an `(n_bands, n_wannier, n_kpoints)` array, or the
`(up, dn)` pair for a [`ProductLayout`](@ref). Intermediates that the gradient
pullback needs again---the `X`/`Y` blocks and assembled `U`---are stashed in
`ws`, so `pullback_gradient!` does not recompute them.
"""
function assemble_gauge! end

"""
    pullback_gradient!(g, layout, model, ws) -> g

Apply the adjoint derivative of gauge assembly to the canonical gradient
`dΩ/dU*`, which the objective left in `ws.GU`, and write the result into the
layout-native gradient buffer `g`. A compact layout omits coordinates fixed by
the frozen subspace, so only active derivatives are stored.
"""
function pullback_gradient! end

# ---- ULayout --------------------------------------------------------------

initial_parameters(::ULayout, model) = copy(model.gauges)

# x ≡ U, so assembly is the identity — no copy, the objective reads `x` directly.
assemble_gauge!(::ULayout, x::AbstractArray{<:Complex, 3}, model, ws) = x

pullback_gradient!(g::AbstractArray{<:Complex, 3}, ::ULayout, model, ws) =
    copyto!(g, ws.GU)

# ---- XYLayout -------------------------------------------------------------

struct _XYBlock
    frozen::Vector{Int}
    nonfrozen::Vector{Int}
    x_range::UnitRange{Int}
    y_range::UnitRange{Int}
    frozen_first::Bool
end

"""Static ranges and row maps for compact `(X,Y)` storage at every kpoint."""
struct _XYStructure
    nbands::Int
    nwann::Int
    blocks::Vector{_XYBlock}
    nparameters::Int
end

function _xy_structure(frozen::AbstractMatrix{Bool}, nwann::Int)
    nbands, nkpts = size(frozen)
    blocks = _XYBlock[]
    offset = 0
    for ik in 1:nkpts
        frozen_rows = findall(view(frozen, :, ik))
        nonfrozen_rows = findall(!, view(frozen, :, ik))
        nfrozen = length(frozen_rows)
        nfrozen <= nwann || error(
            "kpoint $ik has $nfrozen frozen bands but only $nwann Wannier functions"
        )
        x_range = (offset + 1):(offset + nwann^2)
        offset += nwann^2
        ny = length(nonfrozen_rows) * (nwann - nfrozen)
        y_range = (offset + 1):(offset + ny)
        offset += ny
        frozen_first = frozen_rows == 1:nfrozen
        push!(blocks, _XYBlock(frozen_rows, nonfrozen_rows, x_range, y_range, frozen_first))
    end
    return _XYStructure(nbands, nwann, blocks, offset)
end

function _initialize_compact_y!(Y::AbstractArray, xy::_XYStructure)
    fill!(Y, zero(eltype(Y)))
    for (ik, block) in enumerate(xy.blocks)
        for (column, row) in enumerate(block.frozen)
            Y[row, column, ik] = one(eltype(Y))
        end
    end
    return Y
end

# Factor a canonical gauge directly into the compact optimizer coordinates.
# Frozen rows supply the fixed columns of Y; only the complementary Stiefel
# block is stored, so no full Y array is materialized during initialization.
function _initial_xy_parameters(
        U::AbstractArray{T, 3}, frozen::AbstractMatrix{Bool}, xy::_XYStructure
    ) where {T <: Complex}
    nkpts = length(xy.blocks)
    size(U) == (xy.nbands, xy.nwann, nkpts) || throw(
        DimensionMismatch(
            "gauge has size $(size(U)); expected $((xy.nbands, xy.nwann, nkpts))"
        )
    )
    size(frozen) == (xy.nbands, nkpts) || throw(
        DimensionMismatch(
            "frozen mask has size $(size(frozen)); expected $((xy.nbands, nkpts))"
        )
    )

    x = zeros(T, xy.nparameters)
    nwann = xy.nwann
    @inbounds for (ik, block) in enumerate(xy.blocks)
        nfrozen = length(block.frozen)
        nactive = nwann - nfrozen
        Af = orthonormalize_frozen(view(U, :, :, ik), view(frozen, :, ik))
        X = reshape(view(x, block.x_range), nwann, nwann)

        nfrozen > 0 && copyto!(view(X, 1:nfrozen, :), view(Af, block.frozen, :))
        if nactive > 0
            Ar = view(Af, block.nonfrozen, :)
            P = Ar * Ar'
            vectors = eigen(Hermitian((P + P') / 2)).vectors
            Yactive = reshape(
                view(x, block.y_range), length(block.nonfrozen), nactive
            )
            copyto!(Yactive, view(vectors, :, (size(vectors, 2) - nactive + 1):size(vectors, 2)))
            mul!(view(X, (nfrozen + 1):nwann, :), Yactive', Ar)
        end
        X .= lowdin_orthonormalize(X)
    end
    return x
end

function _unpack_xy!(X::AbstractArray, Y::AbstractArray, x::AbstractVector, xy::_XYStructure)
    length(x) == xy.nparameters || throw(
        DimensionMismatch(
            "compact XY vector has length $(length(x)); expected $(xy.nparameters)"
        )
    )
    nwann = xy.nwann
    for (ik, block) in enumerate(xy.blocks)
        copyto!(view(X, :, :, ik), reshape(view(x, block.x_range), nwann, nwann))
        nfrozen = length(block.frozen)
        Yactive = reshape(
            view(x, block.y_range), length(block.nonfrozen), nwann - nfrozen
        )
        if block.frozen_first
            copyto!(
                view(Y, (nfrozen + 1):xy.nbands, (nfrozen + 1):nwann, ik),
                Yactive,
            )
        else
            for column in axes(Yactive, 2), (row_index, row) in enumerate(block.nonfrozen)
                Y[row, nfrozen + column, ik] = Yactive[row_index, column]
            end
        end
    end
    return X, Y
end

function _form_u_compact!(U::AbstractArray, X::AbstractArray, Y::AbstractArray, xy::_XYStructure)
    nwann = xy.nwann
    for (ik, block) in enumerate(xy.blocks)
        nfrozen = length(block.frozen)
        nactive = nwann - nfrozen
        if block.frozen_first
            nfrozen > 0 && copyto!(
                view(U, 1:nfrozen, :, ik), view(X, 1:nfrozen, :, ik)
            )
            if nactive > 0
                mul!(
                    view(U, (nfrozen + 1):xy.nbands, :, ik),
                    view(Y, (nfrozen + 1):xy.nbands, (nfrozen + 1):nwann, ik),
                    view(X, (nfrozen + 1):nwann, :, ik),
                )
            elseif xy.nbands > nfrozen
                fill!(view(U, (nfrozen + 1):xy.nbands, :, ik), zero(eltype(U)))
            end
            continue
        end

        for (column, row) in enumerate(block.frozen), n in 1:nwann
            U[row, n, ik] = X[column, n, ik]
        end
        for row in block.nonfrozen, n in 1:nwann
            value = zero(eltype(U))
            for column in 1:nactive
                value += Y[row, nfrozen + column, ik] * X[nfrozen + column, n, ik]
            end
            U[row, n, ik] = value
        end
    end
    return U
end

function assemble_gauge!(U, X, Y, x::AbstractVector, xy::_XYStructure)
    _unpack_xy!(X, Y, x, xy)
    return _form_u_compact!(U, X, Y, xy)
end

function pullback_gradient!(
        g::AbstractVector, GU::AbstractArray, X::AbstractArray, Y::AbstractArray,
        xy::_XYStructure,
    )
    length(g) == xy.nparameters || throw(
        DimensionMismatch(
            "compact XY gradient has length $(length(g)); expected $(xy.nparameters)"
        )
    )
    nwann = xy.nwann
    for (ik, block) in enumerate(xy.blocks)
        nfrozen = length(block.frozen)
        nactive = nwann - nfrozen
        GX = reshape(view(g, block.x_range), nwann, nwann)
        GY = reshape(view(g, block.y_range), length(block.nonfrozen), nactive)

        if block.frozen_first
            nfrozen > 0 && copyto!(view(GX, 1:nfrozen, :), view(GU, 1:nfrozen, :, ik))
            if nactive > 0
                mul!(
                    view(GX, (nfrozen + 1):nwann, :),
                    view(Y, (nfrozen + 1):xy.nbands, (nfrozen + 1):nwann, ik)',
                    view(GU, (nfrozen + 1):xy.nbands, :, ik),
                )
                mul!(
                    GY,
                    view(GU, (nfrozen + 1):xy.nbands, :, ik),
                    view(X, (nfrozen + 1):nwann, :, ik)',
                )
            end
            continue
        end

        for row in 1:nfrozen, n in 1:nwann
            GX[row, n] = GU[block.frozen[row], n, ik]
        end
        for row in 1:nactive, n in 1:nwann
            value = zero(eltype(g))
            for band in block.nonfrozen
                value += conj(Y[band, nfrozen + row, ik]) * GU[band, n, ik]
            end
            GX[nfrozen + row, n] = value
        end
        for column in 1:nactive, (row_index, band) in enumerate(block.nonfrozen)
            value = zero(eltype(g))
            for n in 1:nwann
                value += GU[band, n, ik] * conj(X[nfrozen + column, n, ik])
            end
            GY[row_index, column] = value
        end
    end
    return g
end

function initial_parameters(::XYLayout, model)
    xy = _xy_structure(model.frozen_bands, n_wannier(model))
    return _initial_xy_parameters(model.gauges, model.frozen_bands, xy)
end

function assemble_gauge!(::XYLayout, x::AbstractVector, model, ws)
    return assemble_gauge!(ws.U, ws.X, ws.Y, x, ws.xy)
end

pullback_gradient!(g::AbstractVector, ::XYLayout, model, ws) =
    pullback_gradient!(g, ws.GU, ws.X, ws.Y, ws.xy)

# ---- ProductLayout -------------------------------------------------------

function initial_parameters(::ProductLayout{XYLayout, XYLayout}, model)
    up = _xy_structure(model.up.frozen_bands, n_wannier(model.up))
    dn = _xy_structure(model.dn.frozen_bands, n_wannier(model.dn))
    return vcat(
        _initial_xy_parameters(model.up.gauges, model.up.frozen_bands, up),
        _initial_xy_parameters(model.dn.gauges, model.dn.frozen_bands, dn),
    )
end

function assemble_gauge!(::ProductLayout{XYLayout, XYLayout}, x::AbstractVector, model, ws)
    nup = ws.up.xy.nparameters
    return (
        assemble_gauge!(ws.up.U, ws.up.X, ws.up.Y, view(x, 1:nup), ws.up.xy),
        assemble_gauge!(
            ws.dn.U, ws.dn.X, ws.dn.Y, view(x, (nup + 1):length(x)), ws.dn.xy
        ),
    )
end

function pullback_gradient!(
        g::AbstractVector, ::ProductLayout{XYLayout, XYLayout}, model, ws
    )
    nup = ws.up.xy.nparameters
    pullback_gradient!(
        view(g, 1:nup), ws.up.GU, ws.up.X, ws.up.Y, ws.up.xy
    )
    pullback_gradient!(
        view(g, (nup + 1):length(g)), ws.dn.GU, ws.dn.X, ws.dn.Y, ws.dn.xy
    )
    return g
end

# ---- WLayout -------------------------------------------------------------

# A single rotation shared by all kpoints; the starting point is the identity.
initial_parameters(::WLayout, model) =
    Matrix{eltype(model.gauges)}(I, n_wannier(model), n_wannier(model))

# The canonical gauge is the same rotation applied at every kpoint: UW_k = U_k W.
function assemble_gauge!(::WLayout, W::AbstractMatrix, model, ws)
    @inbounds for ik in axes(ws.U, 3)
        mul!(view(ws.U, :, :, ik), view(model.gauges, :, :, ik), W)
    end
    return ws.U
end

# Chain rule for UW_k = U_k W collapses the k index: dΩ/dW* = Σ_k U_k† dΩ/dUW_k*.
function pullback_gradient!(g::AbstractMatrix, ::WLayout, model, ws)
    fill!(g, 0)
    @inbounds for ik in axes(ws.GU, 3)
        mul!(g, view(model.gauges, :, :, ik)', view(ws.GU, :, :, ik), true, true)
    end
    return g
end

"""
    finalize_result(layout, x, model)

Convert the final optimizer parameters `x` into the caller-facing result. Gauge
layouts return canonical gauges with no aliasing of workspace buffers;
[`WLayout`](@ref) returns its optimized rotation matrix directly.
"""
function finalize_result end

finalize_result(::ULayout, x, model) = x
finalize_result(::WLayout, x, model) = x

function finalize_result(::XYLayout, x, model)
    xy = _xy_structure(model.frozen_bands, n_wannier(model))
    T = eltype(x)
    X = zeros(T, xy.nwann, xy.nwann, length(xy.blocks))
    Y = zeros(T, xy.nbands, xy.nwann, length(xy.blocks))
    U = similar(Y)
    _initialize_compact_y!(Y, xy)
    return assemble_gauge!(U, X, Y, x, xy)
end

function finalize_result(::ProductLayout{XYLayout, XYLayout}, x, model)
    up = _xy_structure(model.up.frozen_bands, n_wannier(model.up))
    nup = up.nparameters
    return finalize_result(XYLayout(), view(x, 1:nup), model.up),
        finalize_result(XYLayout(), view(x, (nup + 1):length(x)), model.dn)
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

struct _XYManifoldBlock
    x_range::UnitRange{Int}
    y_range::UnitRange{Int}
    nwann::Int
    ynrows::Int
    yncols::Int
end

"""Product of the variably sized Stiefel factors in a compact `XYLayout`."""
struct XYManifold <: Optim.Manifold
    blocks::Vector{_XYManifoldBlock}
end

_shift_range(range::UnitRange, offset::Int) =
    (first(range) + offset):(last(range) + offset)

function XYManifold(structures::_XYStructure...)
    blocks = _XYManifoldBlock[]
    offset = 0
    for xy in structures
        for block in xy.blocks
            nfrozen = length(block.frozen)
            push!(
                blocks,
                _XYManifoldBlock(
                    _shift_range(block.x_range, offset),
                    _shift_range(block.y_range, offset),
                    xy.nwann,
                    length(block.nonfrozen),
                    xy.nwann - nfrozen,
                ),
            )
        end
        offset += xy.nparameters
    end
    return XYManifold(blocks)
end

function Optim.retract!(M::XYManifold, x)
    for block in M.blocks
        X = reshape(view(x, block.x_range), block.nwann, block.nwann)
        X .= lowdin_orthonormalize(X)
        if !isempty(block.y_range)
            Y = reshape(view(x, block.y_range), block.ynrows, block.yncols)
            Y .= lowdin_orthonormalize(Y)
        end
    end
    return x
end

function Optim.project_tangent!(M::XYManifold, g, x)
    for block in M.blocks
        X = reshape(view(x, block.x_range), block.nwann, block.nwann)
        GX = reshape(view(g, block.x_range), block.nwann, block.nwann)
        GX .-= X * ((X' * GX .+ GX' * X) ./ 2)
        if !isempty(block.y_range)
            Y = reshape(view(x, block.y_range), block.ynrows, block.yncols)
            GY = reshape(view(g, block.y_range), block.ynrows, block.yncols)
            GY .-= Y * ((Y' * GY .+ GY' * Y) ./ 2)
        end
    end
    return g
end

function manifold(::XYLayout, model)
    return XYManifold(_xy_structure(model.frozen_bands, n_wannier(model)))
end

manifold(::WLayout, _model) = Optim.Stiefel_SVD()

function manifold(::ProductLayout{XYLayout, XYLayout}, model)
    up = _xy_structure(model.up.frozen_bands, n_wannier(model.up))
    dn = _xy_structure(model.dn.frozen_bands, n_wannier(model.dn))
    return XYManifold(up, dn)
end
