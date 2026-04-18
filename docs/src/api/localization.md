# Localization

The localization API is built around three orthogonal concepts:

- **`Objective`** — the scalar functional to minimize (Marzari-Vanderbilt variance, variance + WF-center penalty, co-optimized spin-polarized spread, …). One concrete subtype per localization variant.
- **`Layout`** — how the gauge parameter is packed on the Stiefel manifold (`UGauge`, `XYGauge`, `ProductLayout`, `WLayout`).
- **`AbstractLocalizationSolver`** — the optimizer backend. `OptimLBFGS` (the default) drives `Optim.jl`; additional backends (e.g. `Manopt.jl`) can be added without touching the objective or layout code.

A `Problem` bundles `(objective, model, layout, workspace)` and is constructed per optimization run. `solve!(problem, solver)` returns the optimized gauge without mutating the `Model`.

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
U_up, U_dn = localize(sm; λs = 1.0)

# Rotation-only refinement (single W matrix)
W = localize(Variance(), model, WLayout())
```

Under the hood every one of these calls expands to

```julia
solve!(Problem(objective, model, layout), OptimLBFGS(; kwargs...))
```

## Migration from the pre-rewrite API

```julia
# before                                   # after
max_localize(model; …)                     localize(model; …)
disentangle(model; …)                      localize(model; …)
disentangle(model, r0, λ; …)               localize(CenteredVariance(r0, λ), model; …)
coopt(sm; λs=1.0, …)                       localize(sm; λs=1.0, …)
constrain_center_coopt(sm, r0, λ; λs, …)   localize(CenteredCoOptVariance(r0, λ, λs), sm; …)
opt_rotate(model; …)                       localize(Variance(), model, WLayout(); …)

get_fg!_disentangle(model)                 fg! = Wannier._make_optim_fg!(Problem(Variance(), model))
get_fg!_maxloc(model)                      fg! = Wannier._make_optim_fg!(Problem(Variance(), model))
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
