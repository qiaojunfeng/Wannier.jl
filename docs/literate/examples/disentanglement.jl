# # Disentanglement of entangled manifold

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will run the disentanglement algorithm on
the silicon valence + conduction bands.
=#

# ## Preparation
# Load the package
using Wannier
using Wannier.Datasets

#=
## Model construction

We load a pre-computed dataset
=#
model = load_dataset("Si2")

#=
!!! tip

    The [`load_dataset`](@ref) internally calls [`read_w90`](@ref) function, which
    will parse the `win` file and set the frozen window for the `Model` according to
    the `dis_froz_min` and `dis_froz_max` parameters in the `win` file.
    However, you can also change these parameters by calling the
    [`set_frozen_win!`](@ref) function.
=#

#=
## Disentanglement and maximal localization

The [`localize`](@ref) function disentangles and maximally localizes the spread
functional, and returns the final gauge matrices `U`,
=#
U = localize(model);

# The initial spreads are
spread(model)

# The final spreads are
spread(model, U)

#=
It seems contradictory that the final spread is larger than the initial spread,
since the frozen window constraint reduces the degrees of freedom for optimization;
Comparing the final spread with the `Initial spread (with states frozen)` output,
it is clear that the minimization decreases the spread to a large extent.

!!! note

    See keyword arguments of [`localize`](@ref) for convergence thresholds.
=#

#=
## Save the new gauge

Again, we can save the new gauge to an `amn` file,
=#
write_amn("si2.dis.amn", U)

#=
Great! Now you have finished the disentanglement tutorial.

As you may have noticed, the workflow is identical to the previous tutorial —
the same [`localize`](@ref) function drives both the isolated and entangled
cases. Internally it picks a [`ULayout`](@ref) layout when the model has no
frozen/disentangled bands and an [`XYLayout`](@ref) layout when it does; the
caller does not have to distinguish.
=#
