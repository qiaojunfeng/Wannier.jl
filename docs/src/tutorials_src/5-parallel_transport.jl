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

- use random Gaussian as initial projections, then run many iterations of maxiaml localization
- use SCDM for the isolated manifold, then run maximal localization

The PTG does not need any initial projections, rather it constructs smooth
Bloch frames sequentially in the reciprocal space, thus getting rid of
the need of initial projections. Moreover, PTG only works in reciprocal space,
thus more efficient than SCDM, which requires wavefunctions in real space.

## Outline

1. construct a [`Model`](@ref), by reading the `win`, `amn`, `mmn`, and `eig` files
2. truncate the `Model` to the valence bands
3. run `Wannier.jl` [`parallel_transport`](@ref) on the new `Model` to construct PTG
4. further rotate the gauge to remove the gauge arbitraness of the first kpoint
5. maxiaml localize to smoothen the gauge
6. interpolate band, and write the real space WF
7. similarly, Wannierize the top valence band

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

To Wannierize the valence bands, we need to remove the conduction bands from the `model`,
i.e., keeping the first 7 bands (the second `1:7` means targeting 7 WFs)
=#
model_val = Wannier.truncate(model, 1:7, 1:7)
#=
In most cases, the initial gauge after truncation is often useless,
we can inspect the initial spread by
=#
omega(model_val)
#=
not particularly bad but this is mostly due to luck.

## Parallel transport

In such case, we run the [`parallel_transport`](@ref) which will construct a smooth gauge.

Before running PTG, let me assure you that indeed there is no initial projection:
we will generate a series of random unitary matrices, and run PTG starting from
those matrices.

However, to make the tutorial reproducible, we will use a fixed random seed.
=#
using Random
Random.seed!(1234)
R = Wannier.rand_unitary(eltype(model_val.A), size(model_val.A)...)
# and assign it to the `model_val`
model_val.A .= R;
# of course now we have totally destroyed the gauge 😈
omega(model_val)
# and we will ask `parallel_transport` to explicitly use our random matrices,
# by setting argument `use_A = true`,
A2, _ = parallel_transport(model_val; use_A=true)
#=
note the function returns a lit bit extra information which is irrelevant to this tutorial,
so we discard it by a `_`.

Then we inspect the new gauge
=#
omega(model_val, A2)

#=
## Fix gauge at first kpoint

However, the [`parallel_transport`](@ref) itself does not fix the gauge of the first kpoint,
which is arbitrary. We can fix it by maximal localizing w.r.t. a single (i.e. independent of kpoints)
`n_wann * n_wann` rotation matrix `W`, by calling [`opt_rotate`](@ref),
=#
model_val.A .= A2;
W = opt_rotate(model_val)
#=
This is convenient since the calculation is cheap, and helps evade local minimum.

Then let's rotate the input gauge by `W`,
=#
A3 = rotate_A(A2, W);
#=
and inspect the new gauge,
=#
omega(model_val, A3)

#=
## Final maximal localization

Often it is helpful to run a final maximal localization that could further
smoothen the gauge,
=#
model_val.A .= A3;
A4 = max_localize(model_val);
# and inspect the final gauge
omega(model_val, A4)

#=
## Band interpolation

Finally, let's have a look at the band interpolation.

Valence + conduction bands,
=#
model.A .= A;
interp_model = Wannier.InterpolationModel(model)
kpi, E = interpolate(interp_model)
# and top valence band,
model_val.A .= A4;
interp_model_val = Wannier.InterpolationModel(model_val)
kpi2, E2 = interpolate(interp_model_val)
# and plot the band structure
using PlotlyJS
P = Wannier.plot_band_diff(kpi, E, E2)
Main.HTMLPlot(P, 500)  # hide

#=
Oh, no! The band interpolation is not working well 😰.

This is because we are using a very coarse kgrid!
Rerun with a finer kgrid and look at the differences.

## Wannierization of top valence band

Similarly, we can Wannierize only the top valence band, with PTG
=#
model_top = Wannier.truncate(model, [7], [1])
#=
note the `[1]` specifies which WF to be kept, but this does not matter
since we won't use the gauge from the previous `model.A`. We use
this `[1]` to specify that we only need 1 WF.

The initial spread is
=#
omega(model_top)
# with PTG,
A_top, _ = parallel_transport(model_top)
omega(model_top, A_top)
#=
In the single-band case, there is no need to run optimal rotation.
But we can still run maximal localization,
=#
model_top.A .= A_top;
A_top2 = max_localize(model_top)
omega(model_top, A_top2)
# and band interpolation
interp_model_top = Wannier.InterpolationModel(model_top)
kpi3, E3 = interpolate(interp_model_top)
# and compare
P = Wannier.plot_band_diff(kpi, E, E3)
Main.HTMLPlot(P, 500)  # hide

#=
Again, rerun with a finer kgrid and see the improvements of interpolation quality.
Also, you can try to plot the WFs, by first truncating the `unk`
=#
Wannier.truncate_unk(CUR_DIR, [7], "$CUR_DIR/truncate")
# the new `unk`s are stored in a subfolder `truncate`,
# now write the realspace WF
write_realspace_wf("$CUR_DIR/wjl", model_top; unkdir="$CUR_DIR/truncate")
# and visualize with your favorite tool.

#=
Now we have automated Wannierization of the valence manifold,
and the top valence band of MoS2, without any initial projection!
=#
