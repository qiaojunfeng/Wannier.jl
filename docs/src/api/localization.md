# Localization

`localize(...)` has two independent dispatch paths that share no supertype:

- **Gradient-based** — pass a concrete [`Objective`](@ref) (`Variance`, `CenteredVariance`, `SpinCoupledVariance`, `CenteredSpinCoupledVariance`). `Objective` is a scalar functional (mathematical); the call bundles `(objective, model, layout, workspace)` into a [`Problem`](@ref) and hands it to [`solve!`](@ref) with an [`AbstractLocalizationSolver`](@ref Wannier.AbstractLocalizationSolver) backend (`OptimLBFGS` by default; additional backends like `Manopt.jl` plug in here). `Layout` (`ULayout`, `XYLayout`, `ProductLayout`, `WLayout`) picks the parameter packing on the Stiefel manifold.
- **Closed-form** — pass a [`ParallelTransport`](@ref). This routes through [`parallel_transport`](@ref) directly; there is no functional, no `Problem`, no solver.

`Objective` and `ParallelTransport` live at different abstraction levels (a scalar functional versus a whole gauge-construction recipe), so there is no common ancestor. `localize` just dispatches on whichever the caller hands in.

See [Architecture](@ref) for how these types fit together and where to plug in when extending them.

## Quick start

```julia
using Wannier

model = load_dataset("Si2")                     # or read_w90 / etc.

# Everyday entry point — single polymorphic driver.
U = localize(model)                             # → Variance + auto-picked layout
U = localize(model, max_iter = 500, g_tol = 1e-8)

# With a WF-center penalty
U = localize(CenteredVariance(r0, λ), model)

# Co-optimization of spin-polarized WFs on a SpinModel
U_up, U_dn = localize(sm; λ_spin = 1.0)

# Rotation-only refinement (single W matrix)
W = localize(Variance(), model, WLayout())

# Closed-form parallel-transport construction (no solver, no scalar functional)
U = localize(ParallelTransport(), model)
U = localize(ParallelTransport(; use_U = true, log_interp = true), model)
```

Objective calls expand to

```julia
solve!(Problem(objective, model, layout), OptimLBFGS(; kwargs...))
```

To swap the solver, construct a `Problem` explicitly and call `solve!` with the backend of choice:

```julia
using Manopt, Manifolds        # activates the WannierManoptExt extension
prob = Problem(Variance(), model)
U    = solve!(prob, ManoptLBFGS(; g_tol = 1e-8, max_iter = 500))
```

`ParallelTransport` calls expand to [`parallel_transport`](@ref) directly.

!!! note "Manopt.jl backend coverage"
    `ManoptLBFGS` currently implements only the `Variance + ULayout` combination (isolated `max_localize`). Other `(Objective, Layout)` pairs — `XYLayout`, `ProductLayout`, `WLayout` — still require `OptimLBFGS`.

## Migration from the pre-rewrite API

```julia
# before                                   # after
max_localize(model; …)                     localize(model; …)
disentangle(model; …)                      localize(model; …)
disentangle(model, r0, λ; …)               localize(CenteredVariance(r0, λ), model; …)
coopt(sm; λ_spin=1.0, …)                       localize(sm; λ_spin=1.0, …)
constrain_center_coopt(sm, r0, λ; λ_spin, …)   localize(CenteredSpinCoupledVariance(r0, λ, λ_spin), sm; …)
opt_rotate(model; …)                       localize(Variance(), model, WLayout(); …)

get_fg!_disentangle(model)                 fg! = Wannier._optimizer_callback(Problem(Variance(), model))
get_fg!_maxloc(model)                      fg! = Wannier._optimizer_callback(Problem(Variance(), model))
```

The legacy symbol-keyed `LocalizationProblem(:disentangle, …)` / `build_fg!` / `_build_fg_*` / `AbstractLocalizationTerm` / `VarianceTerm` / `CenterConstraintTerm` / `omega(terms, …)` surface no longer exists — every driver now routes through `Problem` + `solve!`.

## Contents

```@contents
Pages = ["localization.md"]
Depth = 2
```

## Index

```@index
Pages = ["localization.md"]
```

## Entry point

```@autodocs
Modules = [Wannier]
Pages   = ["localization/localize.jl"]
```

## Gauge method

```@autodocs
Modules = [Wannier]
Pages   = ["localization/method.jl"]
```

## Objective

```@autodocs
Modules = [Wannier]
Pages   = ["localization/objective.jl"]
```

## Solver

```@autodocs
Modules = [Wannier]
Pages   = ["localization/solver.jl"]
```

## Layouts

```@autodocs
Modules = [Wannier]
Pages   = ["common/layouts.jl"]
```

## Disentanglement helpers

```@autodocs
Modules = [Wannier]
Pages   = ["localization/disentangle.jl"]
```

## Co-optimization (SpinModel)

```@autodocs
Modules = [Wannier]
Pages   = ["localization/coopt.jl"]
```

## Gauge transforms

```@autodocs
Modules = [Wannier]
Pages   = ["localization/gauge.jl"]
```

## Symmetry-adapted (SAWF) localization

Symmetry-constrained localization on the irreducible Brillouin zone composes
with the framework through [`SymmetricModel`](@ref Wannier.SymmetricModel)
and the [`SymmetricXYLayout`](@ref Wannier.SymmetricXYLayout) /
[`SchurLayout`](@ref Wannier.SchurLayout) layouts; see
[Architecture](@ref) for the design.

```@autodocs
Modules = [Wannier]
Pages   = ["symmetry/model.jl", "symmetry/localization.jl"]
```

```@docs
Wannier.clean_littlegroup_reps!
```

## Parallel transport

```@autodocs
Modules = [Wannier]
Pages   = [
    "localization/parallel_transport/parallel_transport.jl",
    "localization/parallel_transport/contraction.jl",
]
```

## Splitting the Model

```@autodocs
Modules = [Wannier]
Pages   = ["localization/split.jl"]
```
