# # 5. Parallel transport for top valence band of MoS2

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will Wannierize the top valence band of MoS2,
by constructing the parallel transport gauge (PTG).

Usually, Wannierization of selected group of bands is difficult because
there is no prior knowledge of the initial projections. Possible solutions are

- use random Gaussian as initial projections, then run may iterations of maxiaml localization
- use SCDM for the isolated manifold, then run maximal localization

The PTG does not need any initial projections, rather it constructs smooth
Bloch frames sequentially in the reciprocal space, thus getting rid of
the need of initial projections. Moreover, PTG only works in reciprocal space,
thus more efficient than SCDM, which requires wavefunctions in real space.

## Outline

1. construct a [`Model`](@ref), by reading the `win`, `amn`, `mmn`, and `eig` files
2. truncate the `Model` to only the top valence band
3. run `Wannier.jl` [`parallel_transport`](@ref) on the new `Model` to construct PTG
4. further rotate the gauge to remove the gauge arbitraness of the first kpoint
5. maxiaml localize to smoothen the gauge
6. interpolate band, and write the real space WF

!!! note

    These tutorials assume you have already been familiar with the
    Wannierization workflow using `QE` and `Wannier90`, a good starting
    point could be the tutorials of
    [`Wannier90`](https://github.com/wannier-developers/wannier90).

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`5-parallel_transport.ipynb`](./5-parallel_transport.ipynb)
    - Julia script: [`5-parallel_transport.jl`](./5-parallel_transport.jl)
=#

# ## Preparation
# Load the package
using Wannier

# Path of current tutorial
CUR_DIR = "5-parallel_transport"

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = read_w90("$CUR_DIR/mos2")

#=
## Maximal localization

Although results from this section is not used in the following sections,
it's good to have a general feeling of the system by maximal localizing
the valence + conduction manifolds of MoS2.

Since we are focusing on the isolated manifolds of MoS2,
we only need to run [`max_localize`](@ref) function,
=#
A = max_localize(model);

# The initial spread is
omega(model)

# The final spread is
omega(model, A)

#=
!!! note

    The convergence thresholds is determined by the
    keyword arguments of [`max_localize`](@ref), e.g., `f_tol` for the tolerance on spread,
    and `g_tol` for the tolerance on the norm of spread gradient, etc. You can use stricter thresholds
    to further minimize a bit the spread.
=#

#=
## Truncate bands

To Wannierize only the top valence band, we need to remove the other bands from the `model`
=#
model2 = Wannier.truncate(model, [7])
#=
After truncatation, we have no initial projection anymore, by default it is filled with
identity matrices. We can inspect the initial spread by
=#
omega(model2)
#=
which shows the gauge is far from optimal.

## Parallel transport

In such case, we run the [`parallel_transport`](@ref) which will construct a smooth gauge
=#
A2, _ = parallel_transport(model2)
#=
note the function returns a lit bit extra information which is irrelevant to this tutorial.

Then we inspect the new gauge
=#
omega(model2, A2)

#=
## Rotate the gauge of the first kpoint

However, the [`parallel_transport`](@ref) itself does not fix the gauge of the first kpoint,
which is arbitrary. We can fix it by maximal localizing w.r.t. a `n_wann * n_wann` single matrix
`W` (i.e. independent of kpoints), by calling [`opt_rotate`](@ref),
=#
model2.A .= A2;
W = opt_rotate(model2)
# then rotate the input gauge by `W`,
A3 = rotate_A(A2, W);
#=
and inspect the new gauge,
=#
omega(model2, A3)

#=
## Final maximal localization to smooth the gauge

Often it is helpful to run a final maximal localization that could further
smoothen the gauge,
=#
model2.A .= A3;
A4 = max_localize(model2)
#=
and inspect the final gauge
=#
omega(model2, A4)

#=
## Band interpolation

Finally, let's have a look at the band interpolation.

Valence + conduction bands,
=#
model.A .= A;
interp_model = Wannier.InterpolationModel(model)
kpi, E = interpolate(interp_model)
#=
and top valence band,
=#
model2.A .= A4;
interp_model2 = Wannier.InterpolationModel(model2)
kpi2, E2 = interpolate(interp_model2)
#=
and plot the band structure
=#
using PlotlyJS
P = Wannier.plot_band_diff(kpi, E, E2)
Main.HTMLPlot(P, 500)  # hide

#=
Now we have automated Wannierization of the top valence band of MoS2,
without any initial projection!
=#
