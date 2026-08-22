# # Under the hood: objectives, layouts, and solvers

#=
```@meta
CurrentModule = Wannier
```
=#

#=
Every localization in the previous tutorials was a one-liner:

```julia
U = localize(model)
```

This tutorial opens that box. `localize` is a thin driver over a handful of
composable types, and knowing them lets you swap the functional being
minimized, the parameterization it is minimized over, and the optimizer that
does the minimizing — independently of each other.

The cast:

| Type | Answers |
|---|---|
| [`Model`](@ref Wannier.Model) | *what system?* — lattice, overlaps, gauges, eigenvalues |
| `Objective` | *what is minimized?* — a scalar functional and its gradient |
| `Layout` | *over what parameters?* — how the gauge is packed for the optimizer |
| [`Problem`](@ref Wannier.Problem) | *this run* — objective + model + layout + scratch buffers |
| `AbstractLocalizationSolver` | *by which optimizer?* — backend, tolerances, linesearch |

See the [Architecture](@ref) page for the design rationale behind this split.
=#

# ## Preparation
using Wannier
using Wannier.Datasets

#=
## The Model: data only

A [`Model`](@ref Wannier.Model) is pure input data. It holds no solver state and no scratch
buffers, so the same `Model` is reused across localization runs,
interpolation, and post-processing.
=#
model = load_dataset("Si2_valence")

#=
Its size is read through accessor functions rather than by poking at fields:
=#
(
    n_bands = n_bands(model),
    n_wannier = n_wannier(model),
    n_kpoints = n_kpoints(model),
    n_bvectors = n_bvectors(model),
    n_atoms = n_atoms(model),
    kgrid_size = kgrid_size(model),
)

#=
The matrices themselves are dense, contiguous arrays — not nested vectors of
matrices — which is what keeps the inner kernels fast:
=#
(overlaps = size(model.overlaps), gauges = size(model.gauges), eigenvalues = size(model.eigenvalues))

#=
Here `n_bands == n_wannier`, so this is an *isolated* manifold:
=#
(isisolated = isisolated(model), isentangled = isentangled(model))

#=
## Objectives: what gets minimized

An `Objective` is a scalar functional plus its gradient. The plain
Marzari–Vanderbilt spread is [`Variance`](@ref Wannier.Variance):
=#
objective = Variance()

#=
Other objectives add terms to it:

- [`CenteredVariance`](@ref Wannier.CenteredVariance)`(r0, λ)` — spread plus a penalty pulling WF
  centers towards targets `r0`; see the
  [Constraining Wannier function centers](@ref) tutorial.
- [`SpinCoupledVariance`](@ref Wannier.SpinCoupledVariance)`(λ_spin)` — spread of both spin channels plus an
  ``\uparrow\downarrow`` overlap coupling, on a [`SpinModel`](@ref Wannier.SpinModel).
- [`CenteredSpinCoupledVariance`](@ref Wannier.CenteredSpinCoupledVariance)`(r0, λ, λ_spin)` — both of the above.

Each is a distinct concrete type rather than a runtime sum of terms, so each
can hand-order its own fused value-and-gradient sweep.

## Layouts: what the optimizer actually varies

The gauge `U` is not always the thing handed to the optimizer. A `Layout`
owns the packing. Which one an objective needs depends on the manifold, and
[`default_layout`](@ref Wannier.default_layout) answers that:
=#
default_layout(objective, model)

#=
For an isolated manifold the parameters are the gauge itself
(`ULayout`). For an entangled manifold the disentanglement introduces a second
block, and the parameters become the packed `XY` matrix (`XYLayout`):
=#
model_entangled = load_dataset("Si2")
(isentangled = isentangled(model_entangled), layout = default_layout(objective, model_entangled))

#=
Layouts also own frozen-band masking and the construction of the matrix
manifold that the optimizer moves on — no other part of the code applies a
frozen mask or assembles a `Stiefel` product by hand.

## Problem: bundling one run

A [`Problem`](@ref Wannier.Problem) ties the pieces together. Construction resolves the
layout and allocates the scratch buffers:
=#
prob = Problem(objective, model)

#=
Unlike the `Model`, a `Problem` is ephemeral — build it, solve it, drop it.
Its `workspace` holds the preallocated buffers (`MU`, `UtMU`, the gradient
accumulator, decoded gauges) that are reused across optimizer iterations
instead of being rebuilt each time:
=#
(layout = prob.layout, workspace = nameof(typeof(prob.workspace)))

#=
You can override the layout when you want something other than the default:

```julia
Problem(objective, model, WLayout())   # optimize a single rotation matrix W
```

## Solvers: which optimizer, and how hard

The solver is a separate object carrying the backend choice, tolerances,
iteration limits, linesearch, and history size:
=#
solver = OptimLBFGS(; max_iter = 4, g_tol = 1.0e-6)

# Keep a copy of the starting gauge, so we can see what `solve!` touches:
gauges_before = copy(model.gauges);

# [`solve!`](@ref Wannier.solve!) runs it and returns the optimized gauge:
U = solve!(prob, solver);

# The spread went down,
spread(model, U)

#=
and — worth emphasizing — `solve!` did **not** touch the model you passed in.
The model still holds the gauge it started with:
=#
model.gauges == gauges_before

#=
That is why `Model` is safe to reuse: nothing that optimizes it mutates it.
Assign the result back explicitly when you do want to keep it, with `.=`
since `Model` is immutable:

```julia
model.gauges .= U
```

## Putting it together

`localize` is exactly the composition of the above:

```julia
localize(obj, model; kwargs...)  ==  solve!(Problem(obj, model), OptimLBFGS(; kwargs...))
```

so these two are the same computation:
=#
U_driver = localize(model; max_iter = 4, g_tol = 1.0e-6)
U_manual = solve!(Problem(Variance(), model), OptimLBFGS(; max_iter = 4, g_tol = 1.0e-6))
U_driver ≈ U_manual

#=
Reach for the explicit form when you want to hold a `Problem` across several
solver settings, or to drive a non-default backend.

## Swapping the layout: rotation-only optimization

Passing a layout explicitly changes what is optimized without changing the
functional. `WLayout` optimizes a single rotation matrix ``W`` shared by all
``k``-points, rather than an independent gauge at every ``k``-point:
=#
W = localize(Variance(), model, WLayout(); max_iter = 4)
size(W)

#=
So the result is one small ``n_{\textrm{wann}} \times n_{\textrm{wann}}``
matrix instead of one gauge per ``k``-point. Fold it into a gauge with
`merge_gauge`:
=#
spread(model, merge_gauge(model.gauges, W))

#=
For this dataset the spread barely moves: the projections are already at a
rotation optimum, so the best global `W` is essentially the identity — notice
how small the gradient norm was at iteration 0 above. Rotation-only
optimization earns its keep after a disentanglement, where the ``k``-point
gauges are converged but the overall orbital mixing is not.
=#

#=
## Swapping the solver: alternative backends

[`OptimLBFGS`](@ref Wannier.OptimLBFGS) (Optim.jl) is the default and covers every objective and
layout. A second backend, [`ManoptLBFGS`](@ref Wannier.ManoptLBFGS), is available as a package
extension built on Manopt.jl and Manifolds.jl:

```julia
using Manopt, Manifolds        # activates the extension
U = solve!(Problem(Variance(), model), ManoptLBFGS(; g_tol = 1e-8, max_iter = 500))
```

It currently implements the `Variance` + `ULayout` combination; the other
combinations still need `OptimLBFGS`. Because both backends drive the same
gradient code, the results agree to optimizer tolerance.

## Not everything is an optimization

[`ParallelTransport`](@ref Wannier.ParallelTransport) constructs a smooth gauge in closed form — there
is no scalar functional, no `Problem`, and no solver involved. It rides the
same `localize` entry point as a separate dispatch path:

```julia
U = localize(ParallelTransport(), model)
```

See the
[parallel transport tutorial](@ref "Parallel transport for Wannierization of MoS2 top valence band")
for what it does and when to prefer it.

## Where to go next

- [Architecture](@ref) — how these types relate, and the extension points for
  adding your own objective, layout, or solver backend.
- [Localization](@ref) — the full API reference for everything used here.
=#
