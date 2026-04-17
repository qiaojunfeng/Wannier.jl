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
    return spread(ws, state).Ω
end

function spread(ws::Workspace, U::AbstractArray{<:Complex, 3})
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

# -------------------------------------------------------------------------
# CenteredVariance (commit O): variance + WF-center penalty on a Model
# -------------------------------------------------------------------------

"""
    CenteredVariance(r0, λ)

Marzari-Vanderbilt variance with a per-WF center penalty
`Ωc = λ · Σₙ |r_n − r0[n]|²`. `r0` is a length-`n_wann` vector of target
centers (Cartesian, Å); `λ` is the penalty strength.
"""
struct CenteredVariance{T <: Real} <: Objective
    r0::Vector{Vec3{T}}
    λ::T
end

required_layout(::CenteredVariance, model::Model) = isentangled(model) ? XYGauge() : UGauge()

allocate_workspace(::CenteredVariance, model::Model, ::Layout) = Workspace(model)

function fg!(
        G,
        obj::CenteredVariance,
        state::AbstractArray{<:Complex, 3},
        ws::Workspace,
        kstencil,
        overlaps,
    )
    compute_MU_UtMU!(ws, kstencil, overlaps, state)
    if G !== nothing && G !== false
        # Center penalty is applied by the penalty-aware omega_grad! kernel
        # in a single sweep; the factor is folded into the kernel via
        # `center_penalty(r0, λ)`.
        omega_grad!(center_penalty(obj.r0, obj.λ), G, ws, kstencil, overlaps)
    end
    Ωbase = omega!(ws, kstencil, overlaps)
    Ωc = sum(n -> obj.λ * sum((Ωbase.r[n] - obj.r0[n]) .^ 2), eachindex(obj.r0))
    return Ωbase.Ω + Ωc
end

# -------------------------------------------------------------------------
# CoOptVariance / CenteredCoOptVariance (commit P): SpinModel objectives
# -------------------------------------------------------------------------

"""
    SpinWorkspace{T}(up, dn)

Paired workspace for `SpinModel` objectives. Up/down channels get
independent `Workspace{T}` buffers; the `SpinModel.M` Bloch overlap stays
on the model until the Problem/Workspace refactor (Q/R) relocates it.
"""
struct SpinWorkspace{T}
    up::Workspace{T}
    dn::Workspace{T}
end

"""
    CoOptVariance(λs)

Co-optimization of two spin channels: `Ω = Ωup + Ωdn + λs · Ωupdn` where
`Ωupdn = n_wann − tr(|⟨u↑|u↓⟩|²)` is the ↑↓ overlap penalty (see
`omega_updn`). Operates on a `SpinModel`.
"""
struct CoOptVariance{T <: Real} <: Objective
    λs::T
end

required_layout(::CoOptVariance, ::SpinModel) = ProductLayout(XYGauge(), XYGauge())

function allocate_workspace(::CoOptVariance, model::SpinModel, ::Layout)
    return SpinWorkspace(Workspace(model.up), Workspace(model.dn))
end

"""
    CenteredCoOptVariance(r0, λ, λs)

`CoOptVariance` plus a shared-center penalty applied on both spin
channels (see `CenteredVariance`).
"""
struct CenteredCoOptVariance{T <: Real} <: Objective
    r0::Vector{Vec3{T}}
    λ::T
    λs::T
end

required_layout(::CenteredCoOptVariance, ::SpinModel) = ProductLayout(XYGauge(), XYGauge())

function allocate_workspace(::CenteredCoOptVariance, model::SpinModel, ::Layout)
    return SpinWorkspace(Workspace(model.up), Workspace(model.dn))
end

# -------------------------------------------------------------------------
# Problem (commit Q): solver-agnostic bundle (objective, model, layout, ws)
# -------------------------------------------------------------------------

"""
    Problem(objective, model)
    Problem(objective, model, layout)

Solver-agnostic bundle carrying one `Objective`, its `Model` (or
`SpinModel`), the `Layout` dictating parameter packing, and the
preallocated `Workspace`. Constructed per optimization run; reused across
iterations but discarded once the run returns. No solver options inside —
solver choice/tolerances/linesearch live on an
[`AbstractLocalizationSolver`](@ref) passed separately to `solve!`.
"""
struct Problem{O <: Objective, M, L <: Layout, W}
    objective::O
    model::M
    layout::L
    workspace::W
end

function Problem(objective::Objective, model, layout::Layout = required_layout(objective, model))
    ws = allocate_workspace(objective, model, layout)
    return Problem(objective, model, layout, ws)
end
