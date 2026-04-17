# Wannier.jl Modernization Plan

## Priorities

1. **Math cleanliness** — one canonical gradient convention, one kernel per operation, no duplicated near-copies of spread / objective / gradient code.
2. **API cleanliness** — polymorphic `Objective` interface, ephemeral `Problem`, explicit accessor functions, no magic property forwarding.
3. **Performance** — dense arrays from the start, fused `fg!` everywhere, preallocated workspaces, GPU-ready storage layout.

Backward compatibility is not a constraint. Remove legacy names outright; no aliases, no deprecation shims.

## Non-goals (this refactor)

- GPU kernel implementation — architect the seams, defer the kernels.
- Mixed-backend support (CPU + GPU in one run).
- Transitional adapters around `WannierIO` — change `WannierIO.jl` in lockstep.
- Solver package migration (Optim.jl stays as default; pluggable solver API only).
- Parallel-transport / split.jl redesign — can adopt the Objective interface later.

## Key design decisions (from discussion)

1. **Dense arrays throughout.** All per-k-point quantities become contiguous multidim arrays. No nested `Vector{Vector{Matrix{…}}}` anywhere in `Model`, layout parameterizations, or hot kernels. Change `WannierIO.jl` parsers to emit dense arrays directly — no adapter layer.
2. **No `Crystal` struct.** `kstencil` is reciprocal-space and doesn't belong in a real-space bundle. `Model` stays flat.
3. **`SpinModel` replaces `MagModel`.** Two `Model`s with a constructor-time invariant that lattice / atoms / kstencil match. Kstencil duplication (~100 KB) is noise next to overlaps (~GB).
4. **Remove `Base.getproperty` forwarding.** Accessor functions (`kpoints(m)`, `n_kpoints(m)`, `reciprocal_lattice(m)`, `kgrid_size(m)`, …) are the only way to read those quantities. Same functions dispatch on `Model` and `SpinModel`, giving a uniform API across both.
5. **`Objective` is a polymorphic interface; one dedicated subtype per localization variant.** Three methods: `value`, `gradient!`, `fg!`. Two trait methods for coupling: `required_layout(obj, model)`, `allocate_workspace(obj, model, layout)`. Wannier localization has four variants (variance on `Model`; variance+center on `Model`; variance+spin-coupling on `SpinModel`; variance+center+spin-coupling on `SpinModel`) — encode each as its own concrete `Objective` subtype rather than composing at runtime through a `Tuple`. Each subtype's `fg!` explicitly orders `MU`/`UtMU` computation, making fusion straightforward and sidestepping any per-iteration memoization bookkeeping. Truly-different future objectives (e.g., symmetry) are also their own subtypes. A `SumObjective{T<:Tuple}` escape hatch is deferred until a caller genuinely needs ad-hoc composition.
6. **`Problem` is ephemeral and solver-agnostic.** Problem bundles `(objective, model, layout, workspace)` only — no solver options inside. Constructed per optimization run, discarded after. `Model` is reused across optimization runs, interpolation, post-processing.
7. **Solver is pluggable via `AbstractLocalizationSolver`.** `solve!(prob, solver)` decouples optimization backend from Problem. Start with `OptimLBFGS <: AbstractLocalizationSolver` (migration continuity). Add `ManoptLBFGS` later as a second backend — Manopt.jl is the better long-term home for Stiefel optimization but swapping packages mid-rewrite is unnecessary risk. Skip `Optimization.jl` — it's a meta-wrapper useful for apps that want user-selectable backends, but for a library, direct integration is cleaner.
8. **`Workspace` replaces `Cache`.** Same concept, better-scoped name. Rename, and fix the `cache.G = G` reassignment hack during the rename (size `Workspace.G` as `(n_bands, n_wannier, n_kpoints)` at construction; never reassign the field).
9. **Factor-of-2 gradient convention.** Pick the derivation's `df = 2 Re⟨∇f, dx⟩` convention internally. Apply the ½ once at the solver adapter boundary. Document in one place.

## Architecture

```
┌────────────────────┐     ┌─────────────────────────┐
│ Model / SpinModel  │ ◄── │ Problem                 │
│ (physics input,    │     │   objective             │
│  long-lived)       │     │   model                 │
│                    │     │   layout                │
│ - lattice          │     │   workspace             │
│ - atoms            │     └──────────┬──────────────┘
│ - kstencil         │                │
│ - overlaps (dense) │     ┌──────────▼──────────────┐
│ - gauges (dense)   │     │  solve!(problem)        │
│ - eigenvalues      │     │  → encode → Optim → dec │
│ - frozen_bands     │     └─────────────────────────┘
│ - entangled_bands  │
└────────────────────┘
```

### Model (data-only, flat)

```julia
struct Model{T}
    lattice::Mat3{T}
    atom_positions::Vector{Vec3{T}}
    atom_labels::Vector{String}
    kstencil::KspaceStencil{T}
    overlaps::Array{Complex{T}, 4}      # (n_bands, n_bands, n_bvectors, n_kpoints)
    gauges::Array{Complex{T}, 3}        # (n_bands, n_wannier, n_kpoints)
    eigenvalues::Matrix{T}              # (n_bands, n_kpoints)
    frozen_bands::BitMatrix             # (n_bands, n_kpoints)
    entangled_bands::BitMatrix          # (n_bands, n_kpoints)
end
```

### SpinModel

```julia
struct SpinModel{T}
    up::Model{T}
    dn::Model{T}
    # constructor enforces: up.lattice == dn.lattice, up.atom_positions == dn.atom_positions,
    #                       up.atom_labels == dn.atom_labels, up.kstencil == dn.kstencil
end
```

### Accessor functions (replacing getproperty forwarding)

```julia
# On Model
kpoints(m::Model)        = m.kstencil.kpoints
n_kpoints(m::Model)      = n_kpoints(m.kstencil)
n_bvectors(m::Model)     = n_bvectors(m.kstencil)
reciprocal_lattice(m::Model)  = reciprocal_lattice(m.kstencil)
real_lattice(m::Model)   = m.lattice
kgrid_size(m::Model)     = m.kstencil.kgrid_size
kpb_k(m::Model)          = m.kstencil.kpb_k
kpb_G(m::Model)          = m.kstencil.kpb_G
bweights(m::Model)       = m.kstencil.bweights
n_bands(m::Model)        = size(m.gauges, 1)
n_wannier(m::Model)      = size(m.gauges, 2)
n_atoms(m::Model)        = length(m.atom_positions)

# On SpinModel — delegate to up (invariant guarantees match)
kpoints(m::SpinModel)        = kpoints(m.up)
n_kpoints(m::SpinModel)      = n_kpoints(m.up)
reciprocal_lattice(m::SpinModel)  = reciprocal_lattice(m.up)
# ... same pattern for the rest
```

Direct-field accessors (`m.overlaps`, `m.gauges`, `m.eigenvalues`, `m.frozen_bands`, `m.entangled_bands`) stay as dot-syntax since those are genuine fields.

### Layout (parameter packing)

```julia
abstract type Layout end
struct UGauge  <: Layout end                    # x ≡ Array{Complex{T}, 3}, (nb, nw, nk)
struct XYGauge <: Layout end                    # x ≡ Matrix,               (nw² + nb*nw) × nk
struct ProductLayout{L1<:Layout, L2<:Layout} <: Layout
    first::L1
    second::L2
end
struct WLayout <: Layout end                    # single W matrix (opt_rotate)

# Interface (backend-agnostic; dispatches on concrete AbstractArray types)
encode!(x, layout, U, frozen_bands)  :: Nothing
decode!(U, layout, x)                :: Nothing
pack_gradient!(g, layout, GU, frozen_bands) :: Nothing
manifold(layout, model, solver)      :: solver-specific Manifold
```

Frozen-band masking lives inside layout encode/decode — no per-call-site masking anywhere else. `manifold` is solver-parameterized so `OptimLBFGS` and a future `ManoptLBFGS` can return different concrete manifold objects for the same layout.

### Objective

```julia
abstract type Objective end

# Core interface (each Objective is a peer; no bundling)
value(obj, state, ws)        :: Real
gradient!(G, obj, state, ws) :: Nothing
fg!(G, obj, state, ws)       :: Real   # fused; default falls back to value + gradient!

# Trait methods for problem construction
required_layout(obj, model)              :: Layout
allocate_workspace(obj, model, layout)   :: Workspace
```

Concrete types (rename from current `AbstractLocalizationTerm` hierarchy) — one per localization variant:

```julia
struct Variance                 <: Objective end                                       # max_localize / disentangle on Model
struct CenteredVariance{T}      <: Objective; r0::Vector{Vec3{T}}; λ::T; end           # constrain_center on Model
struct CoOptVariance{T}         <: Objective; λs::T; end                               # coopt on SpinModel
struct CenteredCoOptVariance{T} <: Objective; r0::Vector{Vec3{T}}; λ::T; λs::T; end    # constrain_center + coopt on SpinModel

# Future non-MV objectives (e.g. symmetry) are also their own Objective subtypes.
```

No runtime `Tuple` composition, no `WeightedSum`, no `MVSpread{P}` bundling. Each subtype's `fg!` hand-writes the fused compute order (`MU` → `UtMU` → `omega` → `omega_grad` + penalty/coupling increments) so there is no per-iteration memoization bookkeeping to get wrong. Add a `SumObjective{T<:Tuple} <: Objective` for ad-hoc composition.

`state` is the decoded gauge (`Array{Complex{T}, 3}` view for `UGauge`, or layout-specific). `Workspace` (renamed from `Cache`) holds preallocated `MU`, `UtMU`, gradient buffer, and any centers or coupling scratch the subtype needs.

### Problem

```julia
struct Problem{O<:Objective, M, L<:Layout, W}
    objective::O
    model::M
    layout::L
    workspace::W
end

function Problem(objective::Objective, model)
    layout = required_layout(objective, model)
    ws     = allocate_workspace(objective, model, layout)
    return Problem(objective, model, layout, ws)
end
```

No solver options inside Problem — the type is solver-agnostic.

### Solver (pluggable backend)

```julia
abstract type AbstractLocalizationSolver end

struct OptimLBFGS{LS} <: AbstractLocalizationSolver
    f_tol::Float64
    g_tol::Float64
    max_iter::Int
    history_size::Int
    linesearch::LS
end

solve!(prob::Problem, solver::OptimLBFGS) :: typeof(prob.model.gauges)
```

Solver owns: backend choice (Optim.jl now; Manopt.jl later via `ManoptLBFGS <: AbstractLocalizationSolver`), tolerances, iteration limits, linesearch, history size, manifold type translation.

Convenience wrapper for callers who don't care about solver choice:

```julia
solve!(prob; kwargs...) = solve!(prob, OptimLBFGS(; kwargs...))
```

## Out-of-scope files

Unchanged by this refactor:

- `src/interpolation/*` — no optimization, no layout concerns.
- `src/realspace/*`, `src/Datasets.jl`, `src/tools/*`, `src/plot.jl`, `src/symmetry.jl`.
- `src/io/w90/{tb,hr,spn,nnkp,chk}.jl` — minor signature updates only if accessor change bleeds in.
- `src/localization/parallel_transport/*`, `src/localization/split.jl` — can adopt Objective interface later.
- `CrystalBase.jl` — no changes; still reexported.

---

## Commit ladder

Each commit ends with: full test suite green + reference baselines match within tolerance.

### Phase 0 — Baselines (safety net)

**Audit existing tests first.** `test/` already has good coverage with tight reference values — use it. Spot-check before writing anything new:

- [test/localization/disentangle.jl](test/localization/disentangle.jl) — pinned `Ω`, `ΩI`, `ΩOD`, `ΩD`, `ω`, `r` at `atol=1e-7`; finite-difference gradient check.
- [test/localization/max_localize.jl](test/localization/max_localize.jl) — pinned `Ω`, `ΩI`, `ΩOD`, `Ω̃`; finite-difference gradient check.
- [test/localization/opt_rotate.jl](test/localization/opt_rotate.jl), [test/localization/coopt.jl](test/localization/coopt.jl), [test/localization/constrain_center/](test/localization/constrain_center/) — check coverage; these are the paths most likely to have gaps.
- [test/spread.jl](test/spread.jl), [test/bvector.jl](test/bvector.jl) — primitives; likely fine as-is.

**Commit A** — fill gaps only. If `opt_rotate`, `coopt`, or `constrain_center` lack pinned objective values or finite-difference gradient checks, add them at the same tolerance style as `disentangle.jl` / `max_localize.jl`. Seed any RNGs explicitly. This is the parity gate for every subsequent commit; don't duplicate what already exists.

### Phase 1 — Storage, structure, hygiene

Structural conversion first, then math hygiene on top of the new shape.

**Commit B0 — Remove getproperty forwarding, add accessor functions.** ✅ Done.
- Delete `Base.propertynames(::Model)` and `Base.getproperty(::Model, ::Symbol)` in `src/model.jl`.
- Add `kpoints`, `n_kpoints`, `n_bvectors`, `reciprocal_lattice`, `real_lattice`, `kgrid_size`, `kpb_k`, `kpb_G`, `bweights` as functions.
- Update every call site that currently uses `model.kpoints`, `model.recip_lattice`, `model.kgrid_size`, `model.kpb_k`, `model.kpb_G`, `model.bweights`.
- Net diff is large but mechanical.

**Commit B1 — Dense arrays in Model, coordinated with WannierIO.jl.** ✅ Done.
- Change `Model.overlaps` from `Vector{Vector{Matrix{Complex{T}}}}` to `Array{Complex{T}, 4}` with shape `(n_bands, n_bands, n_bvectors, n_kpoints)`.
- Change `Model.gauges` to `Array{Complex{T}, 3}`, `Model.eigenvalues` to `Matrix{T}`, `Model.frozen_bands` / `entangled_bands` to `BitMatrix`.
- In `WannierIO.jl`: update `src/w90/mmn.jl` parser to emit the 4D overlap array directly; update `src/w90/amn.jl`, `src/w90/eig.jl`, `src/w90/chk.jl` to match. Writer side too.
- In Wannier.jl: update `src/io/w90/*.jl` construction path; update `src/spread.jl` kernels (`compute_MU_UtMU!`, `omega!`, `omega_grad!`) to take slices via `@view overlaps[:, :, ib, ik]`.
- Also: interpolation pipeline migrated to operate on 3D arrays from WannierIO reducers — `simplify` now dispatches on 3D / scalar-matrix / static-array-matrix operator inputs; `fourier` / `invfourier` extended to accept `Array{T,3}`; `eigen!` dense overloads added; `compute_D_matrix` / `projectability` / hamiltonian hessian updated for Matrix eigenvalues + 3D gauges; tests updated to new return shapes.
- Coordinate release: Wannier.jl bumps WannierIO.jl compat lower bound.
- Large blast radius; parity gate is critical here.

**Commit B2 — Replace `MagModel` with `SpinModel`.**
- Rename the type, define it as `{up::Model, dn::Model}`, add constructor invariant checks (`up.lattice == dn.lattice`, `up.atom_positions == dn.atom_positions`, `up.kstencil == dn.kstencil`).
- Delete `MagModel` from `src/localization/coopt.jl` and all direct references.
- Update every file that uses `MagModel` (`localization/coopt.jl`, `localization/constrain_center/coopt.jl`).
- Add accessor dispatches for `SpinModel` (delegate to `up`).

✅ Done. Deviations to follow up later:
- `SpinModel` still carries the Bloch `M` overlap field. The plan target (`up`, `dn` only) requires moving `M` into the coopt `Objective`/`Workspace`, which is deferred until commits P/Q introduce that infrastructure.
- Invariant checks use `isapprox` rather than `==` on `lattice` / `atom_positions` / `kstencil` to tolerate float round-trip through parsers.

**Commit C — Delete broken branches.** ✅ Done.
- Remove `random_gauge` branch in `disentangle.jl` (undefined symbols).
- Other dead code surfaced in review.

**Commit D — Merge `Spread` and `SpreadCenter`.** ✅ Done.
- One `Spread` type with optional center fields (`Nothing` or values). Printing branches on presence.
- Delete `SpreadCenter` and its `# TODO` comment.

**Commit E — Rename `Cache` → `Workspace`, fix `.G` reassignment.** ✅ Done.
- Rename [mutable struct Cache{T}](src/spread.jl#L126) to `Workspace{T}`; update all call sites.
- Fix the `cache.G = G` reassignment pattern in [src/localization/problem.jl](src/localization/problem.jl#L112-L133): size `Workspace.G` as `(n_bands, n_wannier, n_kpoints)` at construction, never reassign the field. Sub-functions that need a local gradient accumulator get a separate buffer, not field mutation.

**Commit F — Factor-of-2 convention.** ✅ Partially done.
- Internal convention: `df = 2 Re⟨∇f, dx⟩` throughout derivation and `omega_grad!`.
- Apply `½` once at the solver adapter boundary (the `fg!` closure passed to Optim).
- Remove the scattered convention notes in `coopt.jl`.

Done so far: collapsed the scattered convention commentary into a single top-of-module note in [src/spread.jl](src/spread.jl) and tightened the `coopt.jl` callout. Deferred: the explicit ½ rescale at the solver-adapter boundary — it requires re-running the parity gate against Optim's LBFGS because the current omission has been absorbed into the line-search behavior, and the explicit rescale is best introduced at the same time as the pluggable solver adapter (commit R).

**Commit G — Collapse `GU_to_GX_GY` and `GU_to_G!`.** ✅ Done.
- Deleted the unused `GU_to_G!` helper in `disentangle.jl`; `GU_to_GX_GY` is the single implementation.

**Commit H — Unify `compute_MU_UtMU!` dispatch.** ✅ Done.
- One dense method over `Array{T,4}` overlaps + `Array{T,3}` gauges; Workspace shim drops the `getfield` work-around now that `Workspace` is a plain immutable struct.

**Commit I — Fuse `SpinModel` objective/gradient.** ✅ Done.
- `_build_fg_mag_disentangle` now produces a single fused `fg!(F, G, XY)` that shares the X/Y decode and up/dn omega paths; `f` and `g!` are thin wrappers returned alongside for backwards-compatible destructuring. `disentangle(::SpinModel, ...)` drives Optim via `Optim.only_fg!(fg!)` instead of independent `f`/`g!` closures.

### Phase 2 — Layout abstraction

**Commit J** — define `Layout`, `UGauge`, `XYGauge`, `ProductLayout`, `WLayout`. Port the following into methods on these layouts: `U_to_X_Y`, `X_Y_to_U`, `XY_to_X_Y`, `X_Y_to_XY`, `GU_to_G!`. One module owns all conversions. ✅ Done (shim layer).
- Types + `encode!`/`decode!`/`pack_gradient!` interface added in [src/common/layouts.jl](src/common/layouts.jl); for now they delegate to the existing `U_to_X_Y` / `X_Y_to_XY` / `GU_to_GX_GY` helpers. Legacy free functions stay callable; objective/problem migration to the layout interface lands with commits Q/R.
- `XYGauge.pack_gradient!` intentionally errors and points callers at `pack_gradient_xy!(g, GU, X, Y, frozen)` since the caller already has decoded X/Y in hand.

**Commit K** — centralize frozen-band masking inside layout encode/decode. Remove per-call-site masking. ✅ Verified.
- Audit confirms masking only fires inside `U_to_X_Y` (XYGauge encode) and `GU_to_GX_GY` / `pack_gradient_xy!` (XYGauge gradient pack) plus the new `UGauge.pack_gradient!`. The remaining `frozen_bands` uses in src are mask *construction* (`get_frozen_bands`, `io/w90/model.jl`, `io/truncate.jl`, `io/w90/chk.jl`) — not per-call-site mask application. The test-only `zero_froz_grad!(G, frozen)` reference-gradient shim stays, as it zeros the NLSolversBase FD gradient that bypassed our layout path. No production call sites apply frozen outside the layout interface.

**Commit L** — move manifold construction (`Stiefel`, `Stiefel_SVD`, `ProductManifold`, `PowerManifold`) into `manifold(layout, model)` methods. Solver code stops constructing manifolds. ✅ Done.
- `manifold(layout, model)` methods live in [src/common/layouts.jl](src/common/layouts.jl); `max_localize`, `disentangle`, `opt_rotate`, `coopt`, and `constrain_center/coopt` call `manifold(UGauge()/XYGauge()/WLayout()/ProductLayout(...), model)` instead of hand-rolling `Optim.ProductManifold(...)` piles.
- Solver-parameterized dispatch (e.g. `manifold(layout, model, ::OptimLBFGS)`) lands with commit R when `AbstractLocalizationSolver` arrives.

Tests (with Commit J):
- Round-trip: `decode(encode(U)) ≈ U`, `encode(decode(x)) ≈ x`, isolated and entangled.
- `pack_gradient!` preserves inner product under finite differences against a known functional.
- `ProductLayout` equals composition of components for a `SpinModel`.
- Frozen-band invariance: decoding preserves frozen rows of `U`.

### Phase 3 — Objective interface

**Commit M** — define `Objective`, `Workspace` contracts (`value`, `gradient!`, `fg!`, `required_layout`, `allocate_workspace`). Empty shell + doctests. ✅ Done.
- New [src/localization/objective.jl](src/localization/objective.jl) introduces `abstract type Objective` and the five contract functions plus a generic `fg!` fallback that routes through `gradient!` + `value`. Concrete subtypes land in N / O / P; the existing `VarianceTerm` / `CenterConstraintTerm` + `:symbol`-keyed `LocalizationProblem` keep running until the migration lands.

**Commit N** — implement `Variance <: Objective`. `required_layout(::Variance, m::Model) = any_entangled(m) ? XYGauge() : UGauge()`. Port existing `omega` / `omega_grad` math into `fg!`. ✅ Done (shell).
- `Variance` subtype, `required_layout` (UGauge/XYGauge branch on `isentangled(model)`), `allocate_workspace(::Variance, model, ::Layout) = Workspace(model)`, and a `fg!` that fuses `compute_MU_UtMU!` + `omega_grad!` + `omega!` live in [src/localization/objective.jl](src/localization/objective.jl). `value` and `gradient!` are intentionally error paths so callers use `fg!` and avoid re-filling MU/UtMU twice; this matches the fusion discipline the plan requires for each Objective subtype.

**Commit O** — implement `CenteredVariance{T} <: Objective` (fields `r0::Vector{Vec3{T}}`, `λ::T`). Same layout trait as `Variance`. `fg!` fuses MV spread with the center penalty in one pass over `MU`/`UtMU`. ✅ Done (shell).
- Concrete `CenteredVariance{T}` lives in [src/localization/objective.jl](src/localization/objective.jl) alongside `Variance`; `required_layout` / `allocate_workspace` mirror the variance case, and `fg!` threads `center_penalty(r0, λ)` into `omega_grad!` so the penalty is folded in during the same MU/UtMU sweep that computes the base spread.

**Commit P** — implement `CoOptVariance{T} <: Objective` (field `λs::T`) and `CenteredCoOptVariance{T} <: Objective` (fields `r0`, `λ`, `λs`) for `SpinModel`. `required_layout(::_, ::SpinModel) → ProductLayout(...)`. Each subtype's `fg!` hand-orders `omega` on up, `omega` on dn, `Ωupdn`, and (for `Centered…`) the center penalty — a single explicit compute sequence per subtype, no tuple traversal at runtime. ✅ Done (shell; `fg!` wiring lands with commits Q/R).
- `CoOptVariance{T}`, `CenteredCoOptVariance{T}`, and a `SpinWorkspace{T}` (paired up/dn `Workspace`) live in [src/localization/objective.jl](src/localization/objective.jl). `required_layout(::_, ::SpinModel) = ProductLayout(XYGauge(), XYGauge())` and `allocate_workspace` returns the paired buffer. The `fg!` implementations that hand-order `omega_up + omega_dn + Ωupdn + penalty` land alongside the Problem/Solver adapter in commits Q/R so that the fused path drives the actual optimizer call.

Finite-difference gradient tests per objective subtype (added with each commit).

### Phase 4 — Problem + Solver adapter

**Commit Q** — solver-agnostic `Problem{O<:Objective, M, L, W}` with fields `(objective, model, layout, workspace)` only. Outer constructor `Problem(objective, model)` dispatches `required_layout(objective, model)` to pick a layout, then `allocate_workspace`. Drop the current `parameterization::Symbol` dispatch and `solver_options` field from [`LocalizationProblem`](src/localization/problem.jl#L12-L17). ✅ Done (shell).
- New `Problem{O,M,L,W}` + outer constructor in [src/localization/objective.jl](src/localization/objective.jl). Legacy `LocalizationProblem(:symbol)` still drives the optimizers; drop lands with commit X after all entry points migrate.

**Commit R** — define `AbstractLocalizationSolver` + concrete `OptimLBFGS`. Implement `solve!(prob, solver::OptimLBFGS)`: builds `Optim.only_fg!` closure over `fg!(G, prob.objective, decoded_state, prob.workspace)` (single objective call — no tuple summation), pulls manifold from `manifold(prob.layout, prob.model, solver)`, applies solver-boundary ½ factor for the gradient convention. Returns the optimized gauge; does not mutate `prob.model`. ✅ Done (Variance branch).
- `AbstractLocalizationSolver` + `OptimLBFGS` live in [src/localization/solver.jl](src/localization/solver.jl) alongside `solve!(prob::Problem{<:Variance, <:Model}, ::OptimLBFGS)` (delegates to the legacy `_build_fg_maxloc` / `_build_fg_disentangle` closures for now to keep a single gradient path under the parity gate). `CenteredVariance` / `CoOptVariance` / `CenteredCoOptVariance` branches land with T–W when the corresponding entry points migrate.
- ½ solver-boundary factor still deferred (see commit F note); introducing it here required also flipping internal omega_grad!, left for Phase 5.

**Commit S** — parity assertion: run one existing path alongside the new `Problem + solve!` path, assert numerical match. Temporary; removed in Phase 5. ✅ Done.
- `test/localization/parity.jl` runs `disentangle(model)` / `max_localize(model)` beside `solve!(Problem(Variance(), model), OptimLBFGS(...))` and asserts gauge equality at `atol = 1e-12`. Because `solve!` currently delegates to `_build_fg_maxloc` / `_build_fg_disentangle`, the trajectories are byte-identical; parity will tighten when the objective's own `fg!` replaces the legacy builders during Phase 5 migration.

### Phase 5 — Migrate entry points

Top-level functions become one-liners over `(objective, model)`:

```julia
max_localize(model; opts...)                       = solve!(Problem(Variance(),                               model), OptimLBFGS(; opts...))
disentangle(model; opts...)                        = solve!(Problem(Variance(),                               model), OptimLBFGS(; opts...))
localize(model, r0, λ; opts...)                    = solve!(Problem(CenteredVariance(r0, λ),                  model), OptimLBFGS(; opts...))
coopt(sm; λs=1.0, opts...)                         = solve!(Problem(CoOptVariance(λs),                        sm),    OptimLBFGS(; opts...))
constrain_center_coopt(sm, r0, λ; λs=1.0, opts...) = solve!(Problem(CenteredCoOptVariance(r0, λ, λs),         sm),    OptimLBFGS(; opts...))
opt_rotate(model; opts...)                         = solve!(Problem(Variance(), model, WLayout()),                    OptimLBFGS(; opts...))
```

Order: hardest case first.

**Commit T** — migrate `disentangle(model)`. ✅ Done.
**Commit U** — migrate `max_localize(model)`. ✅ Done.
**Commit V** — migrate `opt_rotate(model)` (uses `WLayout`; `solve!` branch deepcopies the model, transforms overlaps, and drives `get_fg!_rotate`). ✅ Done.
**Commit W** — migrate `coopt(sm)` on `SpinModel`, keep `disentangle(sm, λ)` as back-compat shim. ✅ Done.
**Commit X** — migrate `constrain_center_coopt(sm, r0, λ; λs)` + `localize(model, r0, λ)` entry points; delete the parity test from Commit S now that all five entry points route through Problem + solve!. ✅ Done.
- Deferred: removal of the `_build_fg_*` / `LocalizationProblem(:symbol)` plumbing. The Objective-native `fg!` for every subtype is a bigger refactor (planned post-Phase 5); solver adapters still delegate to the legacy closures under the hood so the optimization trajectory stays byte-identical.

### Phase 6 — Naming pass

**Commit Y** — rename pass only, one reviewable diff. Decisions to make during commit:
- Verb-first for user-facing functions (`localize`, `disentangle`, `solve!`).
- Noun types (`Problem`, `MVSpread`, `UGauge`, `SpinModel`).
- `omega` vs `spread` — pick one canonical name (lean `spread` user-facing, keep `omega` as internal variable names where it matches literature).
- Confirm accessor function names (`kpoints`, `n_kpoints`, etc.) are final.

### Phase 7 — Performance and alternative backends

**Commit Z** — benchmark harness. Track `omega!`, `omega_grad!`, `fg!`, end-to-end iterations/sec, workspace allocations. Pin to a fixed test system; commit numbers.

**Commit AA** — micro-optimize hot kernels using dense-array-native operations (batched matmul, avoid per-k loops where possible). Only if benchmarks justify.

**Commit BB** — GPU seam audit. `allocate_workspace(obj, model, layout; backend=CPU())` is the single abstraction point. Verify no objective code uses `Array{…}` concretely (use `AbstractArray{…}`). No GPU implementation yet.

**Commit CC (optional)** — add `ManoptLBFGS <: AbstractLocalizationSolver`. Implement `solve!(prob, ::ManoptLBFGS)` using `Manopt.jl` + `Manifolds.jl`. `manifold(layout, model, ::ManoptLBFGS)` returns `Manifolds.Stiefel` variants. This is strictly additive — `OptimLBFGS` keeps working. Benchmark against `OptimLBFGS`; decision to promote `ManoptLBFGS` as default is a later call, not a Phase 7 deliverable.

### Phase 8 — Docs

**Commit CC** — rewrite the localization section of the docs around `Problem` / `Objective` / `Layout`. One-page migration snippet (not a compat guide): "if you were calling `get_fg!_disentangle(model, λ=1.0)`, write `solve!(Problem(MVSpread(...), model))` instead."

---

## Verification strategy

- **Parity gate** after every code-changing commit: objective value tol < 1e-10, gradient 2-norm tol < 1e-8, gauge Frobenius < 1e-6, iteration count drift < 10%.
- **Finite-difference gradient check** per `Objective` concrete type (Phase 3).
- **Layout round-trip tests** per layout (Phase 2).
- **Accessor consistency tests**: `kpoints(m::Model) == kpoints(SpinModel(m, m).up)` and similar (Phase 1, Commit B2).
- **Benchmark numbers** tracked from Phase 7 onward; regression flagged if > 5% slowdown on any hot path.

---

## Risk log

| Risk | Phase | Mitigation |
|------|-------|-----------|
| WannierIO.jl change breaks downstream consumers | 1 (B1) | Coordinate release; bump semver; announce before merging |
| Dense array entangled case with uniform band count assumption | 1 (B1) | Uniform is already the invariant in current code; add explicit assert |
| Dedicated subtypes duplicate MV spread kernel logic | 3 (N–P) | Factor shared math into module-private helpers (`_omega_body!`, `_omega_grad_body!`); each subtype's `fg!` calls the helper then adds its own penalty / coupling term |
| `WLayout` for opt_rotate doesn't fit general pattern | 5 (U) | Treat as its own layout type; workspace caches pre-rotated overlaps |
| `SpinModel` workspace doubles memory | 3-5 | `ProductLayout` composes two workspaces + one Ωupdn buffer; acceptable |
| Factor-of-2 fix introduces subtle sign/scale bug | 1 (F) | Parity gate + finite-difference check catches it immediately; one-commit revert |
| Accessor function rename churn | 1 (B0) | Mechanical; single pass; covered by test suite |

---

## Commit cadence summary

A → B0, B1, B2, C–I (9 hygiene/structure commits) → J–L (3 layout) → M–P (4 objective) → Q–S (3 problem + solver adapter + parity gate) → T–Y (6 migration) → Z (naming) → AA–CC (3 perf + optional Manopt backend) → DD (docs).

**~32 commits total, each independently revertable and parity-gated.**

Biggest single commit: **B1** (dense arrays + WannierIO coordination). Fallback split: Wannier.jl internal change first with a thin adapter reading nested → dense at load time; WannierIO.jl change + adapter removal second.

## Current state of the rewrite branch (as of plan update)

The `rewrite` branch has already landed: `LocalizationProblem` struct, `AbstractLocalizationTerm` with `VarianceTerm`/`CenterConstraintTerm`, `Cache`, kernel extraction, some layout-ish helpers, magnetic path migrated into the problem struct. Still pending (this plan addresses):

- `Model` is still nested vectors ([src/model.jl:46-65](src/model.jl#L46-L65)) — Commit B1.
- `Base.getproperty` forwarding still present ([src/model.jl:69-86](src/model.jl#L69-L86)) — Commit B0.
- `MagModel` still present, to become `SpinModel` — Commit B2.
- `LocalizationProblem.parameterization::Symbol` dispatch ([src/localization/problem.jl:15-40](src/localization/problem.jl#L15-L40)) — replaced by Layout types in Commits J and Q.
- `cache.G = G` reassignment ([src/localization/problem.jl:112-133](src/localization/problem.jl#L112-L133)) — Commit E.
- `AbstractLocalizationTerm` → `Objective` rename, unify types — Commits M–P.
- `LocalizationProblem.solver_options` field coupled to Optim — Commits Q, R.
- Magnetic paths still use separate `f` + `g!` ([src/localization/problem.jl:230-286](src/localization/problem.jl#L230-L286)) — Commit I.
