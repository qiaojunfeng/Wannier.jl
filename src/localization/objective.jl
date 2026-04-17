"""
Objective interface (commit M shell).

A concrete `Objective <: Objective` bundles the scalar spread functional and
its gradient for one localization variant. Each subtype implements:

    value(obj, state, ws)          :: Real
    gradient!(G, obj, state, ws)   :: Nothing
    fg!(G, obj, state, ws)         :: Real    # fused; default falls back to value + gradient!

plus two trait methods consumed by `Problem`:

    required_layout(obj, model)              :: Layout
    allocate_workspace(obj, model, layout)   :: Workspace

`state` is the decoded canonical gauge (e.g. `Array{Complex{T}, 3}` for
`UGauge`); the layout is responsible for the encode/decode round-trip. `ws`
is the preallocated scratch carrying `MU`, `UtMU`, and any objective-specific
buffers.

The concrete subtypes — `Variance`, `CenteredVariance`, `CoOptVariance`,
`CenteredCoOptVariance` — land in commits N / O / P. This file only sets up
the contract and a generic `fg!` fallback that routes through
`value` + `gradient!`.
"""
abstract type Objective end

"""
    value(obj, state, ws)

Evaluate the scalar spread functional at `state` using the preallocated
workspace `ws`. Subtypes override this.
"""
function value end

"""
    gradient!(G, obj, state, ws)

Write `dΩ/dU*` (or the objective's native gradient) into `G` in place.
"""
function gradient! end

"""
    fg!(G, obj, state, ws)

Fused value + gradient evaluation. Subtypes that can share `MU`/`UtMU`
work between value and gradient should override; the fallback simply
calls `gradient!` followed by `value`.
"""
fg!(G, obj::Objective, state, ws) = (gradient!(G, obj, state, ws); value(obj, state, ws))

"""
    required_layout(obj, model)

Return the [`Layout`](@ref) that `obj` expects the parameter array to use
for this `model`. For `Variance` / `CenteredVariance` this is `UGauge()`
when the manifold is isolated and `XYGauge()` when entangled.
"""
function required_layout end

"""
    allocate_workspace(obj, model, layout)

Construct the preallocated scratch [`Workspace`](@ref) used during
optimization for `(obj, model, layout)`.
"""
function allocate_workspace end

# -------------------------------------------------------------------------
# Variance (commit N): Marzari-Vanderbilt spread, max_localize / disentangle
# -------------------------------------------------------------------------

"""
    Variance()

Marzari-Vanderbilt variance spread functional. Works on `Model` for both
the isolated (`UGauge`) and entangled (`XYGauge`) cases.
"""
struct Variance <: Objective end

required_layout(::Variance, model::Model) = isentangled(model) ? XYGauge() : UGauge()

allocate_workspace(::Variance, model::Model, ::Layout) = Workspace(model)

function value(::Variance, state::AbstractArray{<:Complex, 3}, ws::Workspace)
    return omega(ws, state).Ω
end

function omega(ws::Workspace, U::AbstractArray{<:Complex, 3})
    # recomputes MU/UtMU against the current U; callers that already filled
    # ws.MU/ws.UtMU can call omega!(ws, ...) directly
    error(
        "Variance.value called without ws.MU populated — use fg! which fuses " *
            "compute_MU_UtMU! with omega!/omega_grad!."
    )
end

function gradient!(
        G::AbstractArray{<:Complex, 3},
        ::Variance,
        state::AbstractArray{<:Complex, 3},
        ws::Workspace,
    )
    error(
        "Variance.gradient! called without ws.MU populated — use fg! which " *
            "fuses compute_MU_UtMU! with omega_grad!."
    )
end

function fg!(
        G,
        ::Variance,
        state::AbstractArray{<:Complex, 3},
        ws::Workspace,
        kstencil,
        overlaps,
    )
    compute_MU_UtMU!(ws, kstencil, overlaps, state)
    Ω = nothing
    if G !== nothing && G !== false
        omega_grad!(G, ws, kstencil, overlaps)
    end
    Ω = omega!(ws, kstencil, overlaps).Ω
    return Ω
end
