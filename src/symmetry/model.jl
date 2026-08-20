export SymmetrizedModel, SymXYLayout, SchurLayout

# -----------------------------------------------------------------------------
# Framework integration of symmetry-constrained (SAWF) localization.
#
# `SymmetrizedModel` bundles everything the IBZ-constrained problem needs; the
# two layouts below pack its IBZ variables; and the `Variance` objective
# methods delegate to the existing Level-1/Level-2 kernels of
# src/symmetry/localization.jl. This file is plumbing only — all numerics
# live in the kernels it reuses.
# -----------------------------------------------------------------------------

"""
    SymmetrizedModel(model, constraint, overlaps_ibz)

Model-side bundle for symmetry-constrained (SAWF) localization.

Wraps a full-mesh [`Model`](@ref) whose stencil uses the *global* b-vector
ordering (see [`globalize_stencil`](@ref)) and whose overlaps were unfolded
from the IBZ, together with the [`SymmetryConstraint`](@ref) tables and the
IBZ overlaps `overlaps_ibz` (the `.immn` data, `n_bands × n_bands ×
n_bvectors × n_kpoints_ibz`, global b ordering). The Schur-adapted bases of
[`schur_basis`](@ref) are built lazily on first use (only the
[`SchurLayout`](@ref) needs them).

Works with [`localize`](@ref) / [`Problem`](@ref) + [`solve!`](@ref) like any
other model; the optimization variables are the gauges at the IBZ kpoints
only, and the returned gauge is the `(U_fbz, U_ibz)` pair of the optimized
covariant gauge expanded to the full mesh plus its IBZ representative.
"""
struct SymmetrizedModel{T <: Real}
    model::Model{T}
    constraint::SymmetryConstraint{T}
    overlaps_ibz::Array{Complex{T}, 4}
    # lazily built Schur-adapted bases; `Ref` so the bundle stays immutable
    schur::Base.RefValue{Union{Nothing, SchurBasis{T}}}
end

function SymmetrizedModel(
        model::Model{T},
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
    return SymmetrizedModel{T}(
        model, constraint, Array{Complex{T}, 4}(overlaps_ibz),
        Base.RefValue{Union{Nothing, SchurBasis{T}}}(nothing),
    )
end

# accessors forward to the inner full-mesh model
n_bands(sm::SymmetrizedModel) = n_bands(sm.model)
n_wannier(sm::SymmetrizedModel) = n_wannier(sm.model)
n_kpoints(sm::SymmetrizedModel) = n_kpoints(sm.model)
n_bvectors(sm::SymmetrizedModel) = n_bvectors(sm.model)
n_atoms(sm::SymmetrizedModel) = n_atoms(sm.model)
CrystalBase.real_lattice(sm::SymmetrizedModel) = real_lattice(sm.model)
CrystalBase.reciprocal_lattice(sm::SymmetrizedModel) = reciprocal_lattice(sm.model)
isentangled(sm::SymmetrizedModel) = isentangled(sm.model)
isisolated(sm::SymmetrizedModel) = isisolated(sm.model)

"""Number of IBZ kpoints (the kpoints carrying optimization variables)."""
n_kpoints_ibz(sm::SymmetrizedModel) = sm.constraint.nk_ibz

"""Frozen-band mask restricted to the IBZ kpoints."""
frozen_bands_ibz(sm::SymmetrizedModel) = sm.model.frozen_bands[:, sm.constraint.ibz2fbz]

"""
    schur_basis(sm::SymmetrizedModel)

The lazily built [`SchurBasis`](@ref) of the bundle (built on first call,
cached afterwards).
"""
function schur_basis(sm::SymmetrizedModel)
    sb = sm.schur[]
    sb === nothing || return sb
    sb = schur_basis(sm.constraint, frozen_bands_ibz(sm))
    sm.schur[] = sb
    return sb
end

function Base.show(io::IO, ::MIME"text/plain", sm::SymmetrizedModel)
    println(io, "SymmetrizedModel:")
    println(io, "  n_bands        =  ", n_bands(sm))
    println(io, "  n_wannier      =  ", n_wannier(sm))
    println(io, "  n_kpoints      =  ", n_kpoints(sm), " (FBZ)")
    print(io, "  n_kpoints_ibz  =  ", n_kpoints_ibz(sm))
    return nothing
end

# starting IBZ gauge: the model's IBZ slices projected onto the covariant subspace
function _initial_U_ibz(sm::SymmetrizedModel)
    U0 = extract_ibz_gauges(sm.model.gauges, sm.constraint)
    return project_covariant!(U0, sm.constraint)
end

# -----------------------------------------------------------------------------
# SymXYLayout — (X, Y) blocks at the IBZ kpoints
# -----------------------------------------------------------------------------

"""
    SymXYLayout(level = 2)

[`Layout`](@ref) for [`SymmetrizedModel`](@ref): `x` packs the `(X, Y)`
disentanglement blocks at the IBZ kpoints, as a contiguous
`(n_wannier² + n_bands·n_wannier) × n_kpoints_ibz` matrix. Decoding applies
the covariance projector [`project_covariant!`](@ref); gradient encoding
applies the (self-adjoint) projector followed by the standard `XY` pullback.

`level` selects the objective evaluation the workspace is sized for: `2`
(default) keeps every band-dimension product on the IBZ (the transport
kernels of [`symmetric_fg2!`](@ref)); `1` expands the gauge to the full mesh
each iteration ([`symmetric_fg1!`](@ref)). Identical results, different cost.
"""
struct SymXYLayout <: Layout
    level::Int
    function SymXYLayout(level::Integer = 2)
        level in (1, 2) || error("level must be 1 or 2")
        return new(Int(level))
    end
end

function initial_x(::SymXYLayout, sm::SymmetrizedModel)
    X0, Y0 = U_to_X_Y(_initial_U_ibz(sm), frozen_bands_ibz(sm))
    return X_Y_to_XY(X0, Y0)
end

# Both symmetric workspaces expose the same `U_ibz`/`G_ibz`/`X_ibz`/`Y_ibz`/
# `frozen_ibz` fields, so one decode/encode pair serves Level 1 and Level 2.
function decode!(::SymXYLayout, x::AbstractMatrix, sm::SymmetrizedModel, ws)
    XY_to_X_Y!(ws.X_ibz, ws.Y_ibz, x)
    X_Y_to_U!(ws.U_ibz, ws.X_ibz, ws.Y_ibz)
    return project_covariant!(ws.U_ibz, sm.constraint)
end

function encode_gradient!(g::AbstractMatrix, ::SymXYLayout, sm::SymmetrizedModel, ws)
    project_covariant!(ws.G_ibz, sm.constraint)
    return encode_gradient_xy!(g, ws.G_ibz, ws.X_ibz, ws.Y_ibz, ws.frozen_ibz)
end

function decode(::SymXYLayout, x, sm::SymmetrizedModel)
    sc = sm.constraint
    X, Y = XY_to_X_Y(x, sc.nbands, sc.nwann)
    U_ibz = project_covariant!(X_Y_to_U(X, Y), sc)
    return expand_gauges(U_ibz, sc), U_ibz
end

function manifold(::SymXYLayout, sm::SymmetrizedModel)
    nw, nb = n_wannier(sm), n_bands(sm)
    per_k = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw)
    )
    return Optim.PowerManifold(per_k, (nw^2 + nb * nw,), (n_kpoints_ibz(sm),))
end

# -----------------------------------------------------------------------------
# SchurLayout — flat real Schur-block parameter vector
# -----------------------------------------------------------------------------

"""
    SchurLayout()

[`Layout`](@ref) for [`SymmetrizedModel`](@ref): `x` is the flat REAL vector
of per-(IBZ kpoint, irrep class) Schur block parameters (see
[`schur_basis`](@ref)). Fewer parameters than [`SymXYLayout`](@ref), no
projector calls, and exact little-group covariance (including anti-unitary
elements) by construction; optimized on the [`SchurManifold`](@ref).
"""
struct SchurLayout <: Layout end

"""
Workspace for the [`SchurLayout`](@ref): the Level-2 scratch plus a copy of
the current parameter vector, stashed by `decode!` because the Schur gradient
chain ([`schur_encode_gradient!`](@ref)) needs the parameters again.
"""
struct SchurWorkspace{T}
    inner::SymmetricWorkspace2{T}
    x::Vector{T}
end

initial_x(::SchurLayout, sm::SymmetrizedModel) =
    schur_initial_x(_initial_U_ibz(sm), schur_basis(sm))

function decode!(::SchurLayout, x::AbstractVector, sm::SymmetrizedModel, ws::SchurWorkspace)
    copyto!(ws.x, x)
    return schur_decode!(ws.inner.U_ibz, x, schur_basis(sm))
end

encode_gradient!(g::AbstractVector, ::SchurLayout, sm::SymmetrizedModel, ws::SchurWorkspace) =
    schur_encode_gradient!(g, ws.inner.G_ibz, ws.x, schur_basis(sm))

function decode(::SchurLayout, x, sm::SymmetrizedModel{T}) where {T}
    sc = sm.constraint
    U_ibz = zeros(Complex{T}, sc.nbands, sc.nwann, sc.nk_ibz)
    schur_decode!(U_ibz, x, schur_basis(sm))
    return expand_gauges(U_ibz, sc), U_ibz
end

manifold(::SchurLayout, sm::SymmetrizedModel) = SchurManifold(schur_basis(sm))

# -----------------------------------------------------------------------------
# Variance objective on a SymmetrizedModel
# -----------------------------------------------------------------------------

# the canonical gradient of the constrained problem lives at the IBZ kpoints
_canonical_gradient(ws::SymmetricWorkspace) = ws.G_ibz
_canonical_gradient(ws::SymmetricWorkspace2) = ws.G_ibz
_canonical_gradient(ws::SchurWorkspace) = ws.inner.G_ibz

default_layout(::Union{Variance, CenteredVariance}, ::SymmetrizedModel) = SymXYLayout()

function allocate_workspace(
        ::Union{Variance, CenteredVariance}, sm::SymmetrizedModel,
        layout::SymXYLayout; backend = CPU(),
    )
    layout.level == 1 && return SymmetricWorkspace(sm.model, sm.constraint)
    return SymmetricWorkspace2(sm.model.eigenvalues, sm.model.frozen_bands, sm.constraint)
end

function allocate_workspace(
        ::Union{Variance, CenteredVariance}, sm::SymmetrizedModel{T},
        ::SchurLayout; backend = CPU(),
    ) where {T}
    inner = SymmetricWorkspace2(sm.model.eigenvalues, sm.model.frozen_bands, sm.constraint)
    return SchurWorkspace(inner, zeros(T, schur_basis(sm).nx))
end

# `U_ibz` is the covariant IBZ gauge the layout decoded (normally aliasing the
# workspace buffer the kernels read); `G_ibz` receives the *unprojected*
# canonical gradient dΩ/dU*(ki) — projecting is the layout's job.
function fg!(
        F, G_ibz, ::Variance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetrizedModel,
        ws::SymmetricWorkspace2,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    return _fg2_core!(F, G_ibz, sm.overlaps_ibz, sm.constraint, ws)
end

fg!(
    F, G_ibz, obj::Variance,
    U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetrizedModel, ws::SchurWorkspace,
) = fg!(F, G_ibz, obj, U_ibz, sm, ws.inner)

function fg!(
        F, G_ibz, ::Variance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetrizedModel,
        ws::SymmetricWorkspace,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    return _fg1_core!(F, G_ibz, sm.model, sm.constraint, ws)
end

# -----------------------------------------------------------------------------
# CenteredVariance objective on a SymmetrizedModel
# -----------------------------------------------------------------------------

# The center penalty rides along inside the penalty-aware gradient kernels
# (the same hook as `omega_grad!` in src/spread.jl); its value term is the
# same `Ωc = λ Σₙ |rₙ − r0ₙ|²` as `omega_center`, computed from the WF
# centers the core leaves behind.
_omega_center_value(r, r0, λ) = λ * sum(n -> sum(abs2, r[n] - r0[n]), eachindex(r0))

function fg!(
        F, G_ibz, obj::CenteredVariance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetrizedModel,
        ws::SymmetricWorkspace2,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    pen = center_penalty(obj.r0, obj.λ)
    Ω = _fg2_core!(F, G_ibz, pen, sm.overlaps_ibz, sm.constraint, ws)
    F === nothing && return nothing
    return Ω + _omega_center_value(ws.r, obj.r0, obj.λ)
end

fg!(
    F, G_ibz, obj::CenteredVariance,
    U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetrizedModel, ws::SchurWorkspace,
) = fg!(F, G_ibz, obj, U_ibz, sm, ws.inner)

function fg!(
        F, G_ibz, obj::CenteredVariance,
        U_ibz::AbstractArray{<:Complex, 3}, sm::SymmetrizedModel,
        ws::SymmetricWorkspace,
    )
    U_ibz === ws.U_ibz || copyto!(ws.U_ibz, U_ibz)
    pen = center_penalty(obj.r0, obj.λ)
    Ω = _fg1_core!(F, G_ibz, pen, sm.model, sm.constraint, ws)
    F === nothing && return nothing
    return Ω + _omega_center_value(ws.full.r, obj.r0, obj.λ)
end

"""
    localize(sm::SymmetrizedModel; kwargs...)

Symmetry-constrained (SAWF) localization: minimize the MV spread over
covariant gauges parameterized at the IBZ kpoints only, via [`Variance`](@ref)
+ [`SymXYLayout`](@ref) (Level 2). Returns `(U_fbz, U_ibz)`. `kwargs` forward
to [`OptimLBFGS`](@ref). Pass a layout explicitly for the other variants,
e.g. `localize(Variance(), sm, SchurLayout())` or
`localize(Variance(), sm, SymXYLayout(1))`.
"""
localize(sm::SymmetrizedModel; kwargs...) = localize(Variance(), sm; kwargs...)
