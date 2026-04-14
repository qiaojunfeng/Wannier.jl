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

1. **Dense arrays throughout.** All per-k-point quantities become contiguous multidim arrays. No nested `Vector{Vector{Matrix{…}}}` anywhere in `Model` or hot kernels. Change `WannierIO.jl` parsers to emit dense arrays directly — no adapter layer.
2. **No `Crystal` struct.** `kstencil` is reciprocal-space and doesn't belong in a real-space bundle. `Model` stays flat.
3. **`SpinModel` replaces `MagModel`.** Two `Model`s with a constructor-time invariant that lattice / atoms / kstencil match. Kstencil duplication (~100 KB) is noise next to overlaps (~GB).
4. **Remove `Base.getproperty` forwarding.** Accessor functions (`kpoints(m)`, `n_kpoints(m)`, `reciprocal_lattice(m)`, `kgrid_size(m)`, …) are the only way to read those quantities. Same functions dispatch on `Model` and `SpinModel`, giving a uniform API across both.
5. **`Objective` is a polymorphic interface, not a sum-of-`Term`s.** Three methods: `value`, `gradient!`, `fg!`. Two trait methods for coupling: `required_layout(obj, model)`, `allocate_workspace(obj, model, layout)`. MV-family penalty composition stays internal to `MVSpread` (preserves `MU` / `UtMU` fusion). `WeightedSum` exists only as an escape hatch for future unrelated objectives.
6. **`Problem` is ephemeral; `Model` is long-lived data.** Problem bundles `(objective, model, layout, workspace, solver_opts)`. Constructed per optimization run, discarded after. Model is reused across optimization runs, interpolation, post-processing.
7. **Factor-of-2 gradient convention.** Pick the derivation's `df = 2 Re⟨∇f, dx⟩` convention internally. Apply the ½ once at the solver adapter boundary. Document in one place.

## Architecture

```
┌────────────────────┐     ┌─────────────────────────┐
│ Model / SpinModel  │ ◄── │ Problem                 │
│ (physics input,    │     │   objective             │
│  long-lived)       │     │   layout                │
│                    │     │   workspace             │
│ - lattice          │     │   solver_opts           │
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
struct UGauge  <: Layout end                    # x ≡ Vector{Matrix}, (nb × nw) per k
struct XYGauge <: Layout end                    # x ≡ Matrix,         (nw² + nb*nw) × nk
struct ProductLayout{L1<:Layout, L2<:Layout} <: Layout
    first::L1
    second::L2
end
struct WLayout <: Layout end                    # single W matrix (opt_rotate)

# Interface
encode(layout, U, frozen_bands)  :: x
decode!(U, layout, x, model)     :: Nothing
pack_gradient!(g, layout, GU)    :: Nothing
manifold(layout, model)          :: Optim.Manifold
```

Frozen-band masking lives inside layout encode/decode — no per-call-site masking anywhere else.

### Objective

```julia
abstract type Objective end

# Core interface
value(obj, state, ws)        :: Real
gradient!(G, obj, state, ws) :: Nothing
fg!(G, obj, state, ws)       :: Real   # fused; default falls back to value + gradient!

# Trait methods for problem construction
required_layout(obj, model)              :: Layout
allocate_workspace(obj, model, layout)   :: Workspace
```

Concrete types:

```julia
struct MVSpread{P}     <: Objective; penalty::P; end   # P = Nothing or a penalty object
struct SpinCoupled{O}  <: Objective; inner::O; λ_s; end
struct WeightedSum{T}  <: Objective; parts::T; weights; end
```

`state` is the decoded gauge (a view of the 3D gauge `Array` slice, or a `Vector{SubArray}` per layout). `Workspace` holds preallocated `MU`, `UtMU`, gradient buffer, centers, and whatever the objective needs.

### Problem

```julia
struct Problem{O<:Objective, M, L<:Layout, W, S}
    objective::O
    model::M
    layout::L
    workspace::W
    solver_opts::S
end

function Problem(obj::Objective, model; solver_opts=default_solver_opts())
    layout = required_layout(obj, model)
    ws     = allocate_workspace(obj, model, layout)
    return Problem(obj, model, layout, ws, solver_opts)
end

solve!(prob::Problem) :: typeof(prob.model.gauges)
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

**Commit A** — lock reference outputs under `test/references/`. For each existing entry point, record: final objective value, gradient 2-norm, final gauge (Frobenius tolerance), iteration count (loose band). Paths: `max_localize`, `disentangle` (isolated + entangled), `opt_rotate`, `coopt` (MagModel), `constrain_center/coopt` (MagModel + center). Seed RNGs explicitly. This is the parity gate for every subsequent commit.

### Phase 1 — Storage, structure, hygiene

Structural conversion first, then math hygiene on top of the new shape.

**Commit B0 — Remove getproperty forwarding, add accessor functions.**
- Delete `Base.propertynames(::Model)` and `Base.getproperty(::Model, ::Symbol)` in `src/model.jl`.
- Add `kpoints`, `n_kpoints`, `n_bvectors`, `reciprocal_lattice`, `real_lattice`, `kgrid_size`, `kpb_k`, `kpb_G`, `bweights` as functions.
- Update every call site that currently uses `model.kpoints`, `model.recip_lattice`, `model.kgrid_size`, `model.kpb_k`, `model.kpb_G`, `model.bweights`.
- Net diff is large but mechanical.

**Commit B1 — Dense arrays in Model, coordinated with WannierIO.jl.**
- Change `Model.overlaps` from `Vector{Vector{Matrix{Complex{T}}}}` to `Array{Complex{T}, 4}` with shape `(n_bands, n_bands, n_bvectors, n_kpoints)`.
- Change `Model.gauges` to `Array{Complex{T}, 3}`, `Model.eigenvalues` to `Matrix{T}`, `Model.frozen_bands` / `entangled_bands` to `BitMatrix`.
- In `WannierIO.jl`: update `src/w90/mmn.jl` parser to emit the 4D overlap array directly; update `src/w90/amn.jl`, `src/w90/eig.jl`, `src/w90/chk.jl` to match. Writer side too.
- In Wannier.jl: update `src/io/w90/*.jl` construction path; update `src/spread.jl` kernels (`compute_MU_UtMU!`, `omega!`, `omega_grad!`) to take slices via `@view overlaps[:, :, ib, ik]`.
- Coordinate release: Wannier.jl bumps WannierIO.jl compat lower bound.
- Large blast radius; parity gate is critical here.

**Commit B2 — Replace `MagModel` with `SpinModel`.**
- Rename the type, define it as `{up::Model, dn::Model}`, add constructor invariant checks (`up.lattice == dn.lattice`, `up.atom_positions == dn.atom_positions`, `up.kstencil == dn.kstencil`).
- Delete `MagModel` from `src/localization/coopt.jl` and all direct references.
- Update every file that uses `MagModel` (`localization/coopt.jl`, `localization/constrain_center/coopt.jl`).
- Add accessor dispatches for `SpinModel` (delegate to `up`).

**Commit C — Delete broken branches.**
- Remove `random_gauge` branch in `disentangle.jl` (undefined symbols).
- Other dead code surfaced in review.

**Commit D — Merge `Spread` and `SpreadCenter`.**
- One `Spread` type with optional center fields (`Nothing` or values). Printing branches on presence.
- Delete `SpreadCenter` and its `# TODO` comment.

**Commit E — Fix `Cache.G` reassignment.**
- Size `Cache.G` as `(n_bands, n_wannier, n_kpoints)` at construction, never reassign the field.
- Call sites that expected a flat gradient matrix get a `reshape` view.

**Commit F — Factor-of-2 convention.**
- Internal convention: `df = 2 Re⟨∇f, dx⟩` throughout derivation and `omega_grad!`.
- Apply `½` once at the solver adapter boundary (the `fg!` closure passed to Optim).
- Remove the scattered convention notes in `coopt.jl`.

**Commit G — Collapse `GU_to_GX_GY` and `GU_to_G!`.**
- One implementation. If the documented "5% speedup I don't know how" reappears, investigate.

**Commit H — Unify `compute_MU_UtMU!` dispatch.**
- One method over dense 3D / 4D array views. No more `Vector{Vector{Matrix}}` path after Commit B1.

**Commit I — Fuse `SpinModel` objective/gradient.**
- Replace separate `f` and `g!` closures in the (now `SpinModel`-based) `get_fg!_disentangle` with a single `fg!`. Pure perf win; no API change yet.

### Phase 2 — Layout abstraction

**Commit J** — define `Layout`, `UGauge`, `XYGauge`, `ProductLayout`, `WLayout`. Port the following into methods on these layouts: `U_to_X_Y`, `X_Y_to_U`, `XY_to_X_Y`, `X_Y_to_XY`, `GU_to_G!`. One module owns all conversions.

**Commit K** — centralize frozen-band masking inside layout encode/decode. Remove per-call-site masking.

**Commit L** — move manifold construction (`Stiefel`, `Stiefel_SVD`, `ProductManifold`, `PowerManifold`) into `manifold(layout, model)` methods. Solver code stops constructing manifolds.

Tests (with Commit J):
- Round-trip: `decode(encode(U)) ≈ U`, `encode(decode(x)) ≈ x`, isolated and entangled.
- `pack_gradient!` preserves inner product under finite differences against a known functional.
- `ProductLayout` equals composition of components for a `SpinModel`.
- Frozen-band invariance: decoding preserves frozen rows of `U`.

### Phase 3 — Objective interface

**Commit M** — define `Objective`, `Workspace` contracts. Empty shell + doctests.

**Commit N** — implement `MVSpread{P}` (P defaults to `Nothing`, can hold a center penalty). Implement `required_layout(::MVSpread, ::Model) = any_entangled(m) ? XYGauge() : UGauge()`, `allocate_workspace`, `fg!`. Port existing `omega` / `omega_grad` math into method body.

**Commit O** — implement `SpinCoupled{O}`. `required_layout(::_, ::SpinModel) → ProductLayout(...)`. Port `Ωupdn` term.

**Commit P** — implement `WeightedSum` (minimal implementation, no call site yet; reserved for future non-MV composition).

Finite-difference gradient tests per objective (added with each commit).

### Phase 4 — Problem + solve!

**Commit Q** — `Problem` struct, outer constructor, `solve!`. Optim.LBFGS only; solver options as a struct (`SolverOpts` with `f_tol`, `g_tol`, `max_iter`, `history_size`, linesearch tag). Solver-boundary ½ factor lives inside `solve!`'s `fg!` closure.

**Commit R** — parity assertion: run one existing path alongside the new `Problem` path, assert parity within tolerance. Temporary; removed in Phase 5.

### Phase 5 — Migrate entry points

Top-level functions become one-liners:

```julia
max_localize(model; opts...)          = solve!(Problem(MVSpread(),                     model; solver_opts=opts))
disentangle(model; opts...)           = solve!(Problem(MVSpread(),                     model; solver_opts=opts))
localize(model, penalty; opts...)     = solve!(Problem(MVSpread(penalty),              model; solver_opts=opts))
coopt(spin_model; λs=1.0, opts...)    = solve!(Problem(SpinCoupled(MVSpread(), λs),    spin_model; solver_opts=opts))
opt_rotate(model; opts...)            = solve!(Problem(MVSpread(),                     model; solver_opts=opts, layout=WLayout()))
```

Order: hardest case first.

**Commit S** — migrate `disentangle`.
**Commit T** — migrate `max_localize`.
**Commit U** — migrate `opt_rotate` (requires `WLayout` with `U_k = U_k^0 · W` decode; workspace caches `U_k^0`).
**Commit V** — migrate `coopt` (`SpinModel`).
**Commit W** — migrate `constrain_center/coopt` (`SpinCoupled(MVSpread(CenterPenalty(...)))`).
**Commit X** — delete all `get_fg!_*` functions, delete `src/localization/constrain_center/coopt.jl` duplication, delete the parity assertion from Commit R.

### Phase 6 — Naming pass

**Commit Y** — rename pass only, one reviewable diff. Decisions to make during commit:
- Verb-first for user-facing functions (`localize`, `disentangle`, `solve!`).
- Noun types (`Problem`, `MVSpread`, `UGauge`, `SpinModel`).
- `omega` vs `spread` — pick one canonical name (lean `spread` user-facing, keep `omega` as internal variable names where it matches literature).
- Confirm accessor function names (`kpoints`, `n_kpoints`, etc.) are final.

### Phase 7 — Performance

**Commit Z** — benchmark harness. Track `omega!`, `omega_grad!`, `fg!`, end-to-end iterations/sec, workspace allocations. Pin to a fixed test system; commit numbers.

**Commit AA** — micro-optimize hot kernels using dense-array-native operations (batched matmul, avoid per-k loops where possible). Only if benchmarks justify.

**Commit BB** — GPU seam audit. `allocate_workspace(obj, model, layout; backend=CPU())` is the single abstraction point. Verify no objective code uses `Array{…}` concretely (use `AbstractArray{…}`). No GPU implementation yet.

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
| `WeightedSum` loses MU fusion if used for MV variants | 3 (P) | Don't use it for MV-family; keep penalty as `MVSpread` field |
| `WLayout` for opt_rotate doesn't fit general pattern | 5 (U) | Treat as its own layout type; workspace caches pre-rotated overlaps |
| `SpinModel` workspace doubles memory | 3-5 | `ProductLayout` composes two workspaces + one Ωupdn buffer; acceptable |
| Factor-of-2 fix introduces subtle sign/scale bug | 1 (F) | Parity gate + finite-difference check catches it immediately; one-commit revert |
| Accessor function rename churn | 1 (B0) | Mechanical; single pass; covered by test suite |

---

## Commit cadence summary

A → B0, B1, B2, C–I (9 cleanup commits) → J–L (3 layout) → M–P (4 objective) → Q–R (2 problem) → S–X (6 migration) → Y (naming) → Z–BB (3 perf) → CC (docs).

**~30 commits total, each independently revertable and parity-gated.**

Biggest single commit: **B1** (dense arrays + WannierIO coordination). Fallback split: Wannier.jl internal change first with a thin adapter reading nested → dense at load time; WannierIO.jl change + adapter removal second.
