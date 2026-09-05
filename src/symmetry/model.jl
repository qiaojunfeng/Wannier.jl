export SymmetricModel, SymmetricXYLayout, SchurLayout
export SymmetricEvaluation, FullMeshEvaluation, IBZFactorizedEvaluation

# -----------------------------------------------------------------------------
# Framework integration of symmetry-constrained (SAWF) localization.
#
# `SymmetricModel` bundles everything the IBZ-constrained problem needs; the
# two layouts below pack its IBZ variables; the evaluation objects select
# between the existing full-mesh and IBZ-factorized kernels of
# src/symmetry/localization.jl. This file is plumbing only — all numerics
# live in the kernels it reuses.
# -----------------------------------------------------------------------------

"""
    SymmetricModel(model, constraint, overlaps_ibz)

Model-side bundle for symmetry-constrained (SAWF) localization.

Wraps a full-mesh [`Model`](@ref) whose stencil uses the *global* b-vector
ordering (see [`globalize_bvector_ordering`](@ref)) and whose overlaps were reconstructed
from the IBZ, together with the [`SymmetryConstraint`](@ref) tables and the
IBZ overlaps `overlaps_ibz` (the `.immn` data, `n_bands × n_bands ×
n_bvectors × n_kpoints_ibz`, global b ordering). The Schur-adapted bases of
[`schur_basis`](@ref) are built lazily on first use (only the
[`SchurLayout`](@ref) needs them).

The wrapper is parametric on the wrapped model type (`SymmetricModel{M}`,
today `M = Model`); see the comment on the struct for the composition plan
and its physics caveat.

Works with [`localize`](@ref) / [`Problem`](@ref) + [`solve`](@ref) like any
other model; the optimization variables are the gauges at the IBZ kpoints
only, and the returned gauge is the `(U_fbz, U_ibz)` pair of the optimized
covariant gauge reconstructed on the full mesh plus its IBZ representative.
"""
# The wrapped-model type `M` is a free parameter: symmetrization is a decorator
# orthogonal to the spin axis, so a future `SymmetricModel{SpinModel}` can
# carry one `SymmetryConstraint` per spin channel and compose the per-channel
# symmetry layouts via `ProductLayout`.
#
# Physics caveat: per-channel constraints are valid only when no (antiunitary)
# symmetry operation couples the spin channels; magnetic systems may need a
# spin-space-group constraint acting on the composite (↑, ↓) gauge instead.
struct SymmetricModel{M, T <: Real}
    model::M
    constraint::SymmetryConstraint{T}
    overlaps_ibz::Array{Complex{T}, 4}
    # lazily built Schur-adapted bases; `Ref` so the bundle stays immutable
    schur::Base.RefValue{Union{Nothing, SchurBasis{T}}}
end

function SymmetricModel(
        model,
        constraint::SymmetryConstraint{T},
        overlaps_ibz::AbstractArray{<:Complex, 4},
    ) where {T}
    n_bands(model) == constraint.nbands && n_wannier(model) == constraint.nwann ||
        error("model and symmetry constraint sizes do not match")
    n_kpoints(model) == constraint.nk_fbz ||
        error("model must live on the full mesh of the symmetry constraint")
    nb = constraint.nbands
    size(overlaps_ibz) == (nb, nb, constraint.nbvecs, constraint.nk_ibz) ||
        error("overlaps_ibz must be n_bands × n_bands × n_bvectors × n_kpoints_ibz")
    return SymmetricModel{typeof(model), T}(
        model, constraint, Array{Complex{T}, 4}(overlaps_ibz),
        Base.RefValue{Union{Nothing, SchurBasis{T}}}(nothing),
    )
end

# accessors forward to the inner full-mesh model; generic in the wrapped
# type — any `M` supporting these accessors works
n_bands(sm::SymmetricModel) = n_bands(sm.model)
n_wannier(sm::SymmetricModel) = n_wannier(sm.model)
n_kpoints(sm::SymmetricModel) = n_kpoints(sm.model)
n_bvectors(sm::SymmetricModel) = n_bvectors(sm.model)
n_atoms(sm::SymmetricModel) = n_atoms(sm.model)
CrystalBase.real_lattice(sm::SymmetricModel) = real_lattice(sm.model)
CrystalBase.reciprocal_lattice(sm::SymmetricModel) = reciprocal_lattice(sm.model)
isentangled(sm::SymmetricModel) = isentangled(sm.model)
isisolated(sm::SymmetricModel) = isisolated(sm.model)

"""Number of IBZ kpoints (the kpoints carrying optimization variables)."""
n_kpoints_ibz(sm::SymmetricModel) = sm.constraint.nk_ibz

"""Frozen-band mask restricted to the IBZ kpoints."""
frozen_bands_ibz(sm::SymmetricModel) = sm.model.frozen_bands[:, sm.constraint.ibz2fbz]

"""
    schur_basis(sm::SymmetricModel)

The lazily built [`SchurBasis`](@ref) of the bundle (built on first call,
cached afterwards).
"""
function schur_basis(sm::SymmetricModel)
    sb = sm.schur[]
    sb === nothing || return sb
    sb = schur_basis(sm.constraint, frozen_bands_ibz(sm))
    sm.schur[] = sb
    return sb
end

function Base.show(io::IO, ::MIME"text/plain", sm::SymmetricModel)
    println(io, "SymmetricModel:")
    println(io, "  n_bands        =  ", n_bands(sm))
    println(io, "  n_wannier      =  ", n_wannier(sm))
    println(io, "  n_kpoints      =  ", n_kpoints(sm), " (FBZ)")
    print(io, "  n_kpoints_ibz  =  ", n_kpoints_ibz(sm))
    return nothing
end

# starting IBZ gauge: the model's IBZ slices projected onto the covariant subspace
function _initial_U_ibz(sm::SymmetricModel)
    U0 = extract_ibz_gauges(sm.model.gauges, sm.constraint)
    return project_covariant!(U0, sm.constraint)
end

# -----------------------------------------------------------------------------
# SymmetricXYLayout — (X, Y) blocks at the IBZ kpoints
# -----------------------------------------------------------------------------

"""
Evaluation strategy for a symmetry-constrained localization problem.

The strategy is independent of the gauge [`Layout`](@ref): it changes how the
objective and canonical gauge gradient are evaluated, without changing the
optimizer variables or the map from those variables to the IBZ gauge.
"""
abstract type SymmetricEvaluation end

"""
    FullMeshEvaluation()

Reference evaluation that reconstructs the gauge on the full Brillouin-zone
mesh, applies the ordinary full-mesh objective, and pulls its canonical gauge
gradient back to the irreducible mesh.
"""
struct FullMeshEvaluation <: SymmetricEvaluation end

"""
    IBZFactorizedEvaluation()

Production evaluation that keeps every band-dimensional product on the
irreducible mesh and evaluates full-mesh star contributions through the
symmetry-reconstruction identities.
"""
struct IBZFactorizedEvaluation <: SymmetricEvaluation end

"""
    SymmetricXYLayout()

[`Layout`](@ref) for [`SymmetricModel`](@ref): `x` packs the `(X, Y)`
disentanglement blocks at the IBZ kpoints as one contiguous vector. Every
`X` is stored in full; each `Y` stores only its active
`(n_bands-n_frozen) × (n_wannier-n_frozen)` block. Gauge assembly applies
the covariance projector [`project_covariant!`](@ref); gradient pullback
applies the (self-adjoint) projector before the compact `XY` chain rule.
"""
struct SymmetricXYLayout <: Layout end

function SymmetricXYLayout(path::Symbol)
    throw(
        ArgumentError(
            "SymmetricXYLayout no longer selects the evaluation strategy; " *
                "use SymmetricXYLayout() with evaluation=FullMeshEvaluation() or " *
                "evaluation=IBZFactorizedEvaluation() (received :$path)",
        )
    )
end

function initial_parameters(::SymmetricXYLayout, sm::SymmetricModel)
    frozen = frozen_bands_ibz(sm)
    xy = _xy_structure(frozen, n_wannier(sm))
    return _initial_xy_parameters(_initial_U_ibz(sm), frozen, xy)
end

# Both symmetric workspaces expose the same `U_ibz`/`G_ibz`/`X_ibz`/`Y_ibz`/
# `frozen_ibz` fields, so one assembly/pullback pair serves both evaluation paths.
function assemble_gauge!(::SymmetricXYLayout, x::AbstractVector, sm::SymmetricModel, ws)
    assemble_gauge!(ws.U_ibz, ws.X_ibz, ws.Y_ibz, x, ws.xy)
    return project_covariant!(ws.U_ibz, sm.constraint)
end

function pullback_gradient!(
        g::AbstractVector, ::SymmetricXYLayout, sm::SymmetricModel, ws
    )
    project_covariant!(ws.G_ibz, sm.constraint)
    return pullback_gradient!(g, ws.G_ibz, ws.X_ibz, ws.Y_ibz, ws.xy)
end

function finalize_result(::SymmetricXYLayout, x, sm::SymmetricModel)
    sc = sm.constraint
    xy = _xy_structure(frozen_bands_ibz(sm), sc.nwann)
    T = eltype(x)
    X = zeros(T, sc.nwann, sc.nwann, sc.nk_ibz)
    Y = zeros(T, sc.nbands, sc.nwann, sc.nk_ibz)
    U_ibz = similar(Y)
    _initialize_compact_y!(Y, xy)
    assemble_gauge!(U_ibz, X, Y, x, xy)
    project_covariant!(U_ibz, sc)
    return reconstruct_gauges(U_ibz, sc), U_ibz
end

function manifold(::SymmetricXYLayout, sm::SymmetricModel)
    return XYManifold(_xy_structure(frozen_bands_ibz(sm), n_wannier(sm)))
end

# -----------------------------------------------------------------------------
# SchurLayout — flat real Schur-block parameter vector
# -----------------------------------------------------------------------------

"""
    SchurLayout()

[`Layout`](@ref) for [`SymmetricModel`](@ref): `x` is the flat REAL vector
of per-(IBZ kpoint, irrep class) Schur block parameters (see
[`schur_basis`](@ref)). Fewer parameters than [`SymmetricXYLayout`](@ref), no
projector calls, and exact little-group covariance (including anti-unitary
elements) by construction; optimized on the [`SchurManifold`](@ref).
Irrep matching is deterministic. Inconsistent representation data or failed
magnetic-corepresentation classification raise an error during setup; no
linear-averaging fallback is used to bypass the gauge constraints.
"""
struct SchurLayout <: Layout end

"""
Workspace for the [`SchurLayout`](@ref): the scratch space selected by the
independent [`SymmetricEvaluation`](@ref), plus a copy of the current parameter
vector. The latter is stashed by `assemble_gauge!` because the Schur-gradient
pullback needs the parameters again.
"""
struct SchurWorkspace{T, W}
    inner::W
    x::Vector{T}
end

initial_parameters(::SchurLayout, sm::SymmetricModel) =
    initial_parameters(_initial_U_ibz(sm), schur_basis(sm))

function assemble_gauge!(
        ::SchurLayout, x::AbstractVector, sm::SymmetricModel, ws::SchurWorkspace
    )
    copyto!(ws.x, x)
    return assemble_gauge!(ws.inner.U_ibz, x, schur_basis(sm))
end

pullback_gradient!(
    g::AbstractVector, ::SchurLayout, sm::SymmetricModel, ws::SchurWorkspace
) = pullback_gradient!(g, ws.inner.G_ibz, ws.x, schur_basis(sm))

function finalize_result(::SchurLayout, x, sm::SymmetricModel{<:Any, T}) where {T}
    sc = sm.constraint
    U_ibz = zeros(Complex{T}, sc.nbands, sc.nwann, sc.nk_ibz)
    assemble_gauge!(U_ibz, x, schur_basis(sm))
    return reconstruct_gauges(U_ibz, sc), U_ibz
end

manifold(::SchurLayout, sm::SymmetricModel) = SchurManifold(schur_basis(sm))

# -----------------------------------------------------------------------------
# Variance objective on a SymmetricModel
# -----------------------------------------------------------------------------

# the canonical gradient of the constrained problem lives at the IBZ kpoints
_canonical_gradient(ws::SymmetricFullMeshWorkspace) = ws.G_ibz
_canonical_gradient(ws::SymmetricTransportWorkspace) = ws.G_ibz
_canonical_gradient(ws::SchurWorkspace) = ws.inner.G_ibz

default_layout(::Union{Variance, CenteredVariance}, ::SymmetricModel) = SymmetricXYLayout()

default_evaluation(
    ::Union{Variance, CenteredVariance}, ::SymmetricModel,
    ::Layout,
) = IBZFactorizedEvaluation()

function allocate_workspace(
        ::Union{Variance, CenteredVariance}, sm::SymmetricModel,
        ::SymmetricXYLayout, ::FullMeshEvaluation; backend = CPU(),
    )
    return SymmetricFullMeshWorkspace(sm.model, sm.constraint)
end

function allocate_workspace(
        ::Union{Variance, CenteredVariance}, sm::SymmetricModel,
        ::SymmetricXYLayout, ::IBZFactorizedEvaluation; backend = CPU(),
    )
    return SymmetricTransportWorkspace(sm.model.eigenvalues, sm.model.frozen_bands, sm.constraint)
end

function allocate_workspace(
        obj::Union{Variance, CenteredVariance}, sm::SymmetricModel{<:Any, T},
        ::SchurLayout, evaluation::SymmetricEvaluation; backend = CPU(),
    ) where {T}
    inner = allocate_workspace(obj, sm, SymmetricXYLayout(), evaluation; backend)
    return SchurWorkspace(inner, zeros(T, schur_basis(sm).nx))
end

# `U_ibz` is the covariant IBZ gauge assembled by the layout (normally aliasing the
# workspace buffer the kernels read); `G_ibz` receives the *unprojected*
# canonical gradient dΩ/dU*(ki) — projecting is the layout's job.
function fg!(
        F, G_ibz, ::Variance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetricModel,
        ws::SymmetricTransportWorkspace,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    return _fg_transport_core!(F, G_ibz, sm.overlaps_ibz, sm.constraint, ws)
end

fg!(
    F, G_ibz, obj::Variance,
    U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetricModel, ws::SchurWorkspace,
) = fg!(F, G_ibz, obj, U_ibz, sm, ws.inner)

function fg!(
        F, G_ibz, ::Variance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetricModel,
        ws::SymmetricFullMeshWorkspace,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    return _fg_fullmesh_core!(F, G_ibz, sm.model, sm.constraint, ws)
end

# -----------------------------------------------------------------------------
# CenteredVariance objective on a SymmetricModel
# -----------------------------------------------------------------------------

# The center penalty rides along inside the penalty-aware gradient kernels
# (the same hook as `omega_grad!` in src/spread.jl); its value term is the
# same `Ωc = λ Σₙ |rₙ − r0ₙ|²` as `omega_center`, computed from the WF
# centers the core leaves behind.
_omega_center_value(r, r0, λ) = λ * sum(n -> sum(abs2, r[n] - r0[n]), eachindex(r0))

function fg!(
        F, G_ibz, obj::CenteredVariance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetricModel,
        ws::SymmetricTransportWorkspace,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    pen = center_penalty(obj.r0, obj.λ)
    Ω = _fg_transport_core!(F, G_ibz, pen, sm.overlaps_ibz, sm.constraint, ws)
    F === nothing && return nothing
    return Ω + _omega_center_value(ws.r, obj.r0, obj.λ)
end

fg!(
    F, G_ibz, obj::CenteredVariance,
    U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetricModel, ws::SchurWorkspace,
) = fg!(F, G_ibz, obj, U_ibz, sm, ws.inner)

function fg!(
        F, G_ibz, obj::CenteredVariance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetricModel,
        ws::SymmetricFullMeshWorkspace,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    pen = center_penalty(obj.r0, obj.λ)
    Ω = _fg_fullmesh_core!(F, G_ibz, pen, sm.model, sm.constraint, ws)
    F === nothing && return nothing
    return Ω + _omega_center_value(ws.full.r, obj.r0, obj.λ)
end

"""
    localize(sm::SymmetricModel; evaluation=IBZFactorizedEvaluation(), kwargs...)

Symmetry-constrained (SAWF) localization: minimize the MV spread over
covariant gauges parameterized at the IBZ kpoints only, via [`Variance`](@ref)
+ [`SymmetricXYLayout`](@ref), evaluated with
[`IBZFactorizedEvaluation`](@ref) by default. Returns `(U_fbz, U_ibz)`.
`kwargs` forward to [`OptimLBFGS`](@ref). The layout and evaluation can be
selected independently, for example
`localize(Variance(), sm, SchurLayout())` or
`localize(Variance(), sm, SymmetricXYLayout(); evaluation=FullMeshEvaluation())`.
"""
function localize(
        sm::SymmetricModel;
        evaluation::SymmetricEvaluation = IBZFactorizedEvaluation(),
        kwargs...,
    )
    return localize(Variance(), sm; evaluation, kwargs...)
end

function localize(
        objective::Union{Variance, CenteredVariance},
        sm::SymmetricModel;
        evaluation::SymmetricEvaluation = IBZFactorizedEvaluation(),
        kwargs...,
    )
    layout = default_layout(objective, sm)
    return localize(objective, sm, layout; evaluation, kwargs...)
end

function localize(
        objective::Union{Variance, CenteredVariance},
        sm::SymmetricModel,
        layout::Layout;
        evaluation::SymmetricEvaluation = default_evaluation(objective, sm, layout),
        warmup::Bool = false,
        kwargs...,
    )
    problem = Problem(objective, sm, layout; evaluation)
    return solve(problem, OptimLBFGS(; kwargs...); warmup).solution
end
