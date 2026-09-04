# Architecture

This page describes how the package is put together: what the core types are,
which one owns which responsibility, and where to plug in when extending it.
For the function-by-function reference, see the
[Localization API](@ref Localization) page.

## The two layers

Wannier.jl separates *physics data* from *optimization machinery*.

```
   data layer (long-lived)              optimization layer (per run)
┌────────────────────────┐        ┌──────────────────────────────────┐
│ Model / SpinModel      │ ◄───── │ Problem                          │
│                        │        │   objective :: Objective         │
│  lattice               │        │   model                          │
│  atom_positions        │        │   layout    :: Layout            │
│  atom_labels           │        │   evaluation                     │
│  kstencil              │        │   workspace :: Workspace         │
│  overlaps              │        └────────────────┬─────────────────┘
│  gauges                │        ┌────────────────▼─────────────────┐
│  eigenvalues           │        │ solve!(problem, solver)          │
│  frozen_bands          │        │   initial parameters → optimize  │
│  entangled_bands       │        │   assemble → fg! → pullback      │
└────────────────────────┘        │   solver :: AbstractLocalization…│
                                  └──────────────────────────────────┘
```

A [`Model`](@ref Wannier.Model) is *data only*. It is constructed once (typically by
`read_w90`) and reused across many operations: localization runs,
interpolation, real-space evaluation, post-processing. It carries no solver
state and no scratch buffers.

A [`Problem`](@ref Wannier.Problem) is *ephemeral*. It bundles everything one optimization run
needs, is consumed by [`solve!`](@ref Wannier.solve!), and is discarded afterwards. It holds
no solver options — those live on the solver object passed separately.

## Data layer

### Model

All per-``k``-point quantities are stored as contiguous dense arrays, not
nested vectors of matrices:

| Field | Type | Shape |
|---|---|---|
| `lattice` | `Mat3{T}` | `3 × 3`, columns are lattice vectors (Å) |
| `atom_positions` | `Vector{Vec3{T}}` | length `n_atoms`, fractional |
| `atom_labels` | `Vector{String}` | length `n_atoms` |
| `kstencil` | `KSpaceStencil{T}` | finite-difference ``\mathbf{b}``-vectors |
| `overlaps` | `Array{Complex{T},4}` | `n_bands × n_bands × n_bvectors × n_kpoints` |
| `gauges` | `Array{Complex{T},3}` | `n_bands × n_wannier × n_kpoints` |
| `eigenvalues` | `Matrix{T}` | `n_bands × n_kpoints` (eV) |
| `frozen_bands` | `BitMatrix` | `n_bands × n_kpoints` |
| `entangled_bands` | `BitMatrix` | `n_bands × n_kpoints` |

Dense storage keeps hot kernels cache-friendly, lets them take slices with
`@view overlaps[:, :, ib, ik]`, and leaves the door open for batched and
device-array backends.

`Model` is deliberately flat — there is no nested `Crystal` bundle, because
`kstencil` is a reciprocal-space object and does not belong in a real-space
container.

### Accessors, not property forwarding

Derived quantities are read through functions, never through magic
`getproperty` forwarding:

```julia
n_bands(model)          n_kpoints(model)      n_bvectors(model)
n_wannier(model)        n_atoms(model)
kpoints(model)          kgrid_size(model)     bweights(model)
kpb_k(model)            kpb_G(model)
real_lattice(model)     reciprocal_lattice(model)
```

Genuine fields (`model.overlaps`, `model.gauges`, `model.eigenvalues`,
`model.frozen_bands`, `model.entangled_bands`) are still reached with dot
syntax. The same accessor functions are defined for `SpinModel`, so code
written against them works unchanged for both.

### InterpolationModel

[`InterpolationModel`](@ref Wannier.InterpolationModel) is the long-lived data
model for pointwise Wannier interpolation. Its primitive operators are a typed
named collection on one common real-space domain; the Hamiltonian is the
`:hamiltonian` entry rather than a special storage type. Construction from a
`Model`, `HrDat`, or `TbDat` absorbs Wigner--Seitz or minimum-distance weights
before evaluation.

The external seam stays small:

```julia
interpolation_model = InterpolationModel(model; real_space = MinimumDistance())
result = interpolate(
    interpolation_model, kpoints, (BandEnergy(), BandVelocity())
)
```

Observable recipes declare primitive operators and derived intermediates. An
internal typed plan unions those requirements, then each k-point batch shares
one Fourier phase block, Hamiltonian evaluation, derivative evaluation, and
Hermitian eigendecomposition. The allocation-reusing `Wannier.interpolate!`
implementation writes into views whose final axis is the current k-point
batch; only the public result arrays scale with the complete request. Adding an
observable therefore extends the recipe interface without adding fields to
`InterpolationModel` or creating another observable-specific interpolator.

When symmetry is supplied, [`WannierSymmetry`](@ref Wannier.WannierSymmetry)
is the shared basis-level object: it owns space-group operations, prescribed
Wannier centers, orbital representations, and orbital cell shifts. Construction
uses it to close and project the real-space coefficient orbits once. The hot
interpolation path is unchanged and needs no per-query group average.
`SymmetryConstraint` contains the same object but separately owns the
IBZ/FBZ, little-group, overlap-transport, and parameterization tables needed by
localization.

### SpinModel

[`SpinModel`](@ref Wannier.SpinModel) is a pair of `Model`s plus the Bloch-basis ``\uparrow\downarrow``
overlap used by the co-optimization coupling term:

```julia
struct SpinModel{T}
    up::Model{T}
    dn::Model{T}
    M::Array{Complex{T},3}   # ⟨u_nk^↑ | u_mk^↓⟩, (n_bands, n_bands, n_kpoints)
end
```

The inner constructor enforces that `up` and `dn` describe the same system —
matching lattice, atom positions, atom labels, and ``k``-space stencil.
Duplicating the stencil across the two channels is negligible next to the
overlap arrays.

## Localization layer

### Objective — what is minimized

An `Objective` is a *scalar functional* plus its gradient. Each localization
variant is its own concrete type rather than a runtime composition of terms,
so every one can hand-order its fused value+gradient sweep:

| Type | Minimizes | Runs on |
|---|---|---|
| [`Variance`](@ref Wannier.Variance) | Marzari–Vanderbilt spread ``\Omega`` | `Model`, `SymmetricModel` |
| [`CenteredVariance`](@ref Wannier.CenteredVariance) | ``\Omega`` + WF-center penalty | `Model`, `SymmetricModel` |
| [`SpinCoupledVariance`](@ref Wannier.SpinCoupledVariance) | ``\Omega_\uparrow + \Omega_\downarrow + \lambda_s \Omega_{\uparrow\downarrow}`` | `SpinModel` |
| [`CenteredSpinCoupledVariance`](@ref Wannier.CenteredSpinCoupledVariance) | co-optimization + center penalty | `SpinModel` |

Each subtype supplies one kernel and three traits:

```julia
fg!(F, GU, obj, U, model, ws)                      # fused value + gradient
default_layout(obj, model)                         # → Layout
default_evaluation(obj, model, layout)             # → evaluation strategy
allocate_workspace(obj, model, layout, evaluation) # → Workspace
```

`fg!` works in **canonical coordinates**: `U` is the gauge array (or the
`(up, dn)` pair for a `SpinModel`), and the gradient `dΩ/dU*` goes into `GU`.
The objective never sees the layout. Both slots follow the Optim.jl convention
— pass `nothing` for `GU` to skip the gradient, `nothing` for `F` to skip the
value.

Fusing matters because the value and the gradient share the expensive ``MU``
and ``U^\dagger M U`` products held in `ws`; computing them in one sweep avoids
repeating that work.

### Layout — how gauges are parameterized

A [`Layout`](@ref Wannier.Layout) owns the map from the optimizer's parameter
container `x` to the canonical gauge array `U`, the corresponding gradient
pullback, and the manifold on which `x` lives.

| Layout | Parameter `x` | Used for |
|---|---|---|
| [`ULayout`](@ref Wannier.ULayout) | `U` itself, `(n_bands, n_wannier, n_kpoints)` | isolated manifold |
| [`XYLayout`](@ref Wannier.XYLayout) | flat vector of full `X` and active `(n_bands-n_frozen) × (n_wannier-n_frozen)` `Y` blocks | entangled manifold |
| [`ProductLayout`](@ref Wannier.ProductLayout) | two layouts side by side | `SpinModel` |
| [`WLayout`](@ref Wannier.WLayout) | a single rotation matrix `W` | rotation-only refinement |
| [`SymmetricXYLayout`](@ref Wannier.SymmetricXYLayout) | the same compact `XY` vector at IBZ kpoints only | `SymmetricModel` |
| [`SchurLayout`](@ref Wannier.SchurLayout) | flat real per-irrep Schur block parameters | `SymmetricModel` |

```julia
initial_parameters(layout, model)          # model.gauges → starting x
assemble_gauge!(layout, x, model, ws)      # x → canonical U; cache intermediates
pullback_gradient!(g, layout, model, ws)   # ws.GU → layout-native gradient
finalize_result(layout, x, model)           # final x → caller-facing result
manifold(layout, model)                     # → the manifold to optimize on
```

Two consequences worth knowing:

- **Frozen-band masking lives here.** Nothing outside the layout applies a
  frozen mask; call sites elsewhere never need to think about it.
- **Manifold construction lives here.** Solvers do not hand-assemble
  `Stiefel` / `ProductManifold` / `PowerManifold` piles; they ask the layout.
- **The final result lives here.** Gauge layouts return canonical gauges that do
  not alias workspace buffers; `WLayout` returns the optimized shared rotation `W`.

The layout is normally picked for you: `default_layout(obj, model)` returns
`ULayout()` for an isolated manifold and `XYLayout()` for an entangled one.
Pass a layout explicitly only when you want something else, e.g. `WLayout()`.

### Workspace — preallocated scratch

A `Workspace` holds the buffers reused across optimizer iterations: `MU`,
`UtMU`, the canonical gradient `GU` (that is `dΩ/dU*`, as distinct from the
layout-native `g`), the assembled `X` / `Y` / `U`, and the WF centers `r`. Buffers are sized once at construction and never reassigned, so
the large per-``k``-point arrays are not rebuilt on every iteration.
`SpinWorkspace` pairs two of them plus the ``\uparrow\downarrow`` overlap used
by the co-optimization objectives.

### Problem and Solver

```julia
Problem(objective, model)                 # layout via default_layout
Problem(objective, model, layout)         # explicit layout
Problem(objective, model, layout; evaluation) # independent evaluation strategy
```

Construction resolves the layout and evaluation, allocates the corresponding
workspace, and stores the five pieces. That is all it does — a `Problem` is
solver-agnostic. (`WLayout` is the
one exception: because a shared rotation only makes sense against a model whose
gauge has been folded into the overlaps, its constructor does that transform on
a copy, leaving your model untouched.)

Because the objective, layout, and evaluation choices compose, one bridge
serves every combination after construction has selected the workspace:

```julia
function _optimizer_callback(prob::Problem)
    obj, model, layout, ws = prob.objective, prob.model, prob.layout, prob.workspace
    return function (F, G, x)
        U = assemble_gauge!(layout, x, model, ws)               # layout
        Ω = fg!(F, G === nothing ? nothing : ws.GU, obj, U, model, ws)   # objective
        G === nothing || pullback_gradient!(G, layout, model, ws) # layout
        return Ω
    end
end
```

and one `solve!` drives it for every objective and layout. There is no
per-combination method to write.

Solvers are pluggable behind `AbstractLocalizationSolver`:

| Solver | Backend | Coverage |
|---|---|---|
| [`OptimLBFGS`](@ref Wannier.OptimLBFGS) | Optim.jl | all objectives and layouts (default) |
| [`ManoptLBFGS`](@ref Wannier.ManoptLBFGS) | Manopt.jl + Manifolds.jl, via package extension | `Variance` + `ULayout` |

The solver owns tolerances, iteration limits, linesearch, and history size:

```julia
U = solve!(Problem(Variance(), model), OptimLBFGS(; g_tol = 1e-8, max_iter = 500))
```

`ManoptLBFGS` lives in a package extension, so the main package does not
depend on Manopt.jl; `using Manopt, Manifolds` activates it.

`solve!` returns the optimized gauge and does **not** mutate `prob.model` —
the model you passed in is still the input you started with.

### The `localize` driver

Everyday use goes through one polymorphic entry point:

```julia
localize(model)                                  # Variance, layout auto-picked
localize(sm; λ_spin = 1.0)                           # SpinCoupledVariance on a SpinModel
localize(CenteredVariance(r0, λ), model)         # explicit objective
localize(Variance(), model, WLayout())           # explicit layout
localize(ParallelTransport(), model)             # closed-form, no solver
```

Objective calls expand to `solve!(Problem(obj, model, layout), OptimLBFGS(; kwargs...))`;
keyword arguments forward to the solver.

### Symmetry-adapted WFs ride the same rails

Symmetry-constrained (SAWF) localization is not a separate driver — it is the
same `Objective` × `Layout` × evaluation × solver composition on a different
model bundle.
A [`SymmetricModel`](@ref Wannier.SymmetricModel) wraps a full-mesh `Model` (global-b stencil,
overlaps reconstructed from the IBZ) together with the
[`SymmetryConstraint`](@ref Wannier.SymmetryConstraint) tables and the IBZ overlaps; the
optimization variables live at the IBZ kpoints only. The evaluation strategy
selects either [`IBZFactorizedEvaluation`](@ref Wannier.IBZFactorizedEvaluation),
which keeps the band-dimensional algebra on the IBZ, or
[`FullMeshEvaluation`](@ref Wannier.FullMeshEvaluation), which reconstructs the
gauge and runs the reference full-mesh kernels. Independently, the layout owns
the parameterization and constraint handling:
[`SymmetricXYLayout`](@ref Wannier.SymmetricXYLayout) assembles through the
covariance projector (and pulls the gradient back through its adjoint), while
[`SchurLayout`](@ref Wannier.SchurLayout) parameterizes the covariant gauges exactly by their
per-irrep Schur blocks — fewer real parameters, no projector calls. Because
these axes are separate, either layout can use either evaluation:

```julia
U_fbz, U_ibz = localize(sm) # SymmetricXYLayout + IBZFactorizedEvaluation
U_fbz, U_ibz = localize(Variance(), sm, SchurLayout())
U_fbz, U_ibz = localize(
    Variance(), sm, SchurLayout(); evaluation = FullMeshEvaluation()
)
U_fbz, U_ibz = solve!(Problem(Variance(), sm), OptimLBFGS(; max_iter = 300))
```

`SymmetricModel` is parametric on the wrapped model (`SymmetricModel{M}`,
today `M = Model`): symmetrization is a decorator orthogonal to the spin
axis. A future `SymmetricModel{SpinModel}` composes through the existing
rails — a [`ProductLayout`](@ref Wannier.ProductLayout) of the symmetry layouts with one
[`SymmetryConstraint`](@ref Wannier.SymmetryConstraint) per spin channel, plus one new transport
identity for the ``\uparrow\downarrow`` coupling overlap. The per-channel
constraints are valid only when no (antiunitary) symmetry operation couples
the spin channels; magnetic systems may instead need a spin-space-group
constraint on the composite gauge.

### Parallel transport is a separate path

[`ParallelTransport`](@ref Wannier.ParallelTransport) is a *gauge-construction recipe*, not a scalar
functional: it builds a smooth gauge in closed form, with no `Problem`, no
`fg!`, and no solver. It participates in `localize` as its own dispatch path
and shares no supertype with `Objective` — the two live at different
abstraction levels, and pretending otherwise would only blur both.

## Gradient convention

The derivation behind `omega_grad!` uses the physics convention

```math
\mathrm{d}f(x) = 2\,\mathrm{Re}\langle \nabla f, \mathrm{d}x \rangle
```

which is why explicit factors of 2 and 4 appear inside the kernel. That is a
statement about the derivation, not a mismatch to correct: those factors are
exactly what make the emitted gradient the true derivative of the spread in the
real-inner-product convention Optim.jl consumes. **No rescaling is applied at
the solver boundary, and none is needed** — the finite-difference gradient
checks in the test suite verify it directly, comparing the analytic gradient
elementwise against finite differences of the same value function.

## Extension points

| You want to… | Do this |
|---|---|
| add a new functional (symmetry, custom penalty, …) | define a new `Objective` subtype with one `fg!` method plus `default_layout`, `default_evaluation`, and `allocate_workspace`. It then works with every compatible layout, evaluation, and solver backend |
| add a new parameterization | define a new `Layout` with `initial_parameters` / `assemble_gauge!` / `pullback_gradient!` / `finalize_result` / `manifold`. It then works with every objective |
| add a symmetry evaluation strategy | define a `SymmetricEvaluation` subtype and its `allocate_workspace` methods; layouts and solvers remain unchanged |
| add a new optimizer backend | define `S <: AbstractLocalizationSolver` and `solve!(::Problem, ::S)` |
| run on a device (GPU) | dispatch `allocate_workspace(obj, model, layout, evaluation; backend)` to return device arrays |

The `backend` keyword on `allocate_workspace` is the single seam where array
storage is chosen; `CPU()` is the only backend implemented today. Kernels are
written against `AbstractArray`, so they do not need to change when another
backend lands.

## Design principles

- **One kernel per operation.** No near-duplicate copies of the spread,
  objective, or gradient code — shared math is factored into helpers that each
  objective calls before adding its own penalty or coupling term.
- **One canonical gradient convention**, documented in one place, converted at
  one boundary.
- **Explicit over magic.** Accessor functions instead of property forwarding;
  concrete objective types instead of runtime term composition; layouts
  instead of symbol-keyed dispatch.
- **Dense from the start.** Contiguous arrays everywhere in the hot path.
