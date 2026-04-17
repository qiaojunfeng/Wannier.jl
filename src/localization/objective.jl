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
