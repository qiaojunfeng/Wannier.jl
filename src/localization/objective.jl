export Variance, CenteredVariance, SpinCoupledVariance, CenteredSpinCoupledVariance
export Problem, default_layout

"""
    CPU()

Backend sentinel. `allocate_workspace(obj, model, layout; backend=CPU())`
is the single abstraction point for swapping array storage; a future GPU
backend (e.g. `CUDA()`) would dispatch `allocate_workspace` to return a
`Workspace` parameterized on device arrays. No GPU implementation yet —
the sentinel only encodes the seam.
"""
struct CPU end

"""
A concrete `Objective` is the scalar functional being minimized, together with
its gradient, for one localization variant.

Each subtype implements one kernel,

    fg!(F, GU, obj, U, model, ws) :: Real

plus two traits consumed when a [`Problem`](@ref) is built:

    default_layout(obj, model)             :: Layout
    allocate_workspace(obj, model, layout) :: Workspace

`fg!` works in **canonical coordinates** — it never sees the layout. Pulling
the gradient back to layout-native parameters is [`pullback_gradient!`](@ref)'s job,
which is what keeps the objective axis and the layout axis independent: a new
objective works with every layout, and a new layout works with every objective.

Finite-difference gradient checks for each subtype live in the test suite:
`test/localization/disentangle.jl` (Variance),
`test/localization/constrain_center/disentangle.jl` (CenteredVariance),
`test/localization/coopt.jl` (SpinCoupledVariance), and
`test/localization/constrain_center/coopt.jl` (CenteredSpinCoupledVariance).
"""
abstract type Objective end

"""
    fg!(F, GU, obj, U, model, ws) -> Ω

Fused value-and-gradient evaluation for `obj` in canonical coordinates.

`U` is the gauge array (`n_bands × n_wannier × n_kpoints`); the gradient
`dΩ/dU*` is written into `GU`. Both slots follow the Optim.jl convention: pass
`nothing` for `GU` to skip the gradient, and `nothing` for `F` to skip the
value, in which case `nothing` is returned.

Value and gradient share the expensive `MU` and `UtMU` products held in `ws`,
which is why they are computed in one call rather than two.
"""
function fg! end

"""
    default_layout(obj, model)

Return the [`Layout`](@ref) that `obj` expects the parameter array to use
for this `model`. For `Variance` / `CenteredVariance` this is `ULayout()`
when the manifold is isolated and `XYLayout()` when entangled.
"""
function default_layout end

"""
    allocate_workspace(obj, model, layout; backend=CPU())

Construct the preallocated scratch [`Workspace`](@ref) used during
optimization for `(obj, model, layout)`. `backend` is the single GPU
seam — a future device backend would dispatch here to return a workspace
holding device arrays. No GPU implementation yet.
"""
function allocate_workspace end

# -------------------------------------------------------------------------
# Variance (commit N): Marzari-Vanderbilt spread, max_localize / disentangle
# -------------------------------------------------------------------------

"""
    Variance()

Marzari-Vanderbilt variance spread functional. Works on `Model` for both
the isolated (`ULayout`) and entangled (`XYLayout`) cases.
"""
struct Variance <: Objective end

default_layout(::Variance, model::Model) = isentangled(model) ? XYLayout() : ULayout()

allocate_workspace(::Variance, model::Model, ::Layout; backend = CPU()) =
    Workspace(model)

function fg!(
        F,
        GU,
        ::Variance,
        U::AbstractArray{<:Complex, 3},
        model::Model,
        ws::Workspace,
    )
    kstencil, overlaps = model.kstencil, model.overlaps
    compute_MU_UtMU!(ws, kstencil, overlaps, U)
    GU === nothing || omega_grad!(GU, ws, kstencil, overlaps)
    F === nothing && return nothing
    return omega!(ws, kstencil, overlaps).Ω
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

default_layout(::CenteredVariance, model::Model) = isentangled(model) ? XYLayout() : ULayout()

allocate_workspace(::CenteredVariance, model::Model, ::Layout; backend = CPU()) =
    Workspace(model)

function fg!(
        F,
        GU,
        obj::CenteredVariance,
        U::AbstractArray{<:Complex, 3},
        model::Model,
        ws::Workspace,
    )
    kstencil, overlaps = model.kstencil, model.overlaps
    compute_MU_UtMU!(ws, kstencil, overlaps, U)
    # The center penalty rides along inside the penalty-aware `omega_grad!`
    # kernel, so it costs no extra sweep over the b-vectors.
    GU === nothing || omega_grad!(center_penalty(obj.r0, obj.λ), GU, ws, kstencil, overlaps)
    F === nothing && return nothing
    return omega_center(omega!(ws, kstencil, overlaps); r0 = obj.r0, λ = obj.λ).Ωt
end

# -------------------------------------------------------------------------
# SpinCoupledVariance / CenteredSpinCoupledVariance (commit P): SpinModel objectives
# -------------------------------------------------------------------------

"""
    SpinWorkspace{T}(up, dn, overlaps_updn)

Paired workspace for `SpinModel` objectives. Up/down channels get
independent `Workspace{T}` buffers; `overlaps_updn` holds the Bloch-basis
``\\uparrow\\downarrow`` overlap that the coupling term needs, copied out of
the `SpinModel` when the workspace is allocated.
"""
struct SpinWorkspace{T}
    up::Workspace{T}
    dn::Workspace{T}
    overlaps_updn::Array{Complex{T}, 3}
end

"""
    SpinCoupledVariance(λ_spin)

Co-optimization of two spin channels: `Ω = Ωup + Ωdn + λ_spin · Ωupdn` where
`Ωupdn = n_wann − tr(|⟨u↑|u↓⟩|²)` is the ↑↓ overlap penalty (see
`omega_updn`). Operates on a `SpinModel`.
"""
struct SpinCoupledVariance{T <: Real} <: Objective
    λ_spin::T
end

default_layout(::SpinCoupledVariance, ::SpinModel) = ProductLayout(XYLayout(), XYLayout())

function allocate_workspace(::SpinCoupledVariance, model::SpinModel, ::Layout; backend = CPU())
    return SpinWorkspace(Workspace(model.up), Workspace(model.dn), Array{eltype(model.overlaps_updn), 3}(model.overlaps_updn))
end

"""
    CenteredSpinCoupledVariance(r0, λ, λ_spin)

`SpinCoupledVariance` plus a shared-center penalty applied on both spin
channels (see `CenteredVariance`).
"""
struct CenteredSpinCoupledVariance{T <: Real} <: Objective
    r0::Vector{Vec3{T}}
    λ::T
    λ_spin::T
end

default_layout(::CenteredSpinCoupledVariance, ::SpinModel) = ProductLayout(XYLayout(), XYLayout())

function allocate_workspace(::CenteredSpinCoupledVariance, model::SpinModel, ::Layout; backend = CPU())
    return SpinWorkspace(Workspace(model.up), Workspace(model.dn), Array{eltype(model.overlaps_updn), 3}(model.overlaps_updn))
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

function Problem(objective::Objective, model, layout::Layout = default_layout(objective, model))
    ws = allocate_workspace(objective, model, layout)
    return Problem(objective, model, layout, ws)
end

# `WLayout` optimizes one rotation shared by all kpoints, which only makes sense
# against a model whose gauge has been folded into the overlaps. Doing that here
# rather than in `solve!` keeps every solver method generic, and leaves the
# caller's model untouched.
function Problem(objective::Objective, model::Model, layout::WLayout)
    nw = n_wannier(model)
    n_bands(model) == nw ||
        error("WLayout needs n_bands == n_wannier; run disentanglement first")
    rotated = deepcopy(model)
    rotated.overlaps .= transform_gauge(rotated.overlaps, rotated.kstencil.kpb_k, rotated.gauges)
    rotated.gauges .= identity_gauge(eltype(rotated.gauges), n_kpoints(rotated), nw)
    return Problem(objective, rotated, layout, allocate_workspace(objective, rotated, layout))
end

function Base.show(io::IO, ::MIME"text/plain", prob::Problem)
    println(io, "Problem:")
    println(io, "  objective  =  ", nameof(typeof(prob.objective)))
    println(
        io, "  model      =  ", nameof(typeof(prob.model)), " (", n_bands(prob.model),
        " bands, ", n_wannier(prob.model), " WFs, ", n_kpoints(prob.model), " kpoints)"
    )
    println(io, "  layout     =  ", prob.layout)
    return print(io, "  workspace  =  ", nameof(typeof(prob.workspace)))
end

# --- SpinModel objectives, in canonical coordinates -----------------------
#
# `U` is the `(up, dn)` gauge pair and `GU` the matching gradient pair. The
# ↑↓ coupling gradient is added here, in canonical coordinates; converting the
# sum once is equivalent to converting each term separately because the
# layout's gradient encoding is linear — and it is one conversion cheaper.

function fg!(F, GU, obj::SpinCoupledVariance, U::Tuple, model::SpinModel, ws::SpinWorkspace)
    Uup, Udn = U
    λ = obj.λ_spin
    compute_MU_UtMU!(ws.up, model.up.kstencil, model.up.overlaps, Uup)
    compute_MU_UtMU!(ws.dn, model.dn.kstencil, model.dn.overlaps, Udn)

    if GU !== nothing
        GUup, GUdn = GU
        omega_grad!(GUup, ws.up, model.up.kstencil, model.up.overlaps)
        omega_grad!(GUdn, ws.dn, model.dn.kstencil, model.dn.overlaps)
        if λ != 0
            Cup, Cdn = omega_updn_grad(ws.overlaps_updn, Uup, Udn)
            GUup .+= λ .* Cup
            GUdn .+= λ .* Cdn
        end
    end

    F === nothing && return nothing
    Ωup = omega!(ws.up, model.up.kstencil, model.up.overlaps).Ω
    Ωdn = omega!(ws.dn, model.dn.kstencil, model.dn.overlaps).Ω
    Ωupdn = λ == 0 ? 0.0 : omega_updn(ws.overlaps_updn, Uup, Udn)
    return Ωup + Ωdn + λ * Ωupdn
end

function fg!(F, GU, obj::CenteredSpinCoupledVariance, U::Tuple, model::SpinModel, ws::SpinWorkspace)
    Uup, Udn = U
    λ = obj.λ_spin
    pen = center_penalty(obj.r0, obj.λ)
    compute_MU_UtMU!(ws.up, model.up.kstencil, model.up.overlaps, Uup)
    compute_MU_UtMU!(ws.dn, model.dn.kstencil, model.dn.overlaps, Udn)

    if GU !== nothing
        GUup, GUdn = GU
        omega_grad!(pen, GUup, ws.up, model.up.kstencil, model.up.overlaps)
        omega_grad!(pen, GUdn, ws.dn, model.dn.kstencil, model.dn.overlaps)
        if λ != 0
            Cup, Cdn = omega_updn_grad(ws.overlaps_updn, Uup, Udn)
            GUup .+= λ .* Cup
            GUdn .+= λ .* Cdn
        end
    end

    F === nothing && return nothing
    Ωup = omega_center(
        omega!(ws.up, model.up.kstencil, model.up.overlaps); r0 = obj.r0, λ = obj.λ
    ).Ωt
    Ωdn = omega_center(
        omega!(ws.dn, model.dn.kstencil, model.dn.overlaps); r0 = obj.r0, λ = obj.λ
    ).Ωt
    Ωupdn = λ == 0 ? 0.0 : omega_updn(ws.overlaps_updn, Uup, Udn)
    return Ωup + Ωdn + λ * Ωupdn
end
