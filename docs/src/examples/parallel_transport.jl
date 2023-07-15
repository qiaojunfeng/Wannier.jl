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

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`tutorial.ipynb`](./tutorial.ipynb)
    - Julia script: [`tutorial.jl`](./tutorial.jl)
=#

# ## Preparation
# Load the package
using Wannier
using WannierPlots

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = read_w90("mos2")

#=
## Maximal localization

Although results from this section is not used in the following sections,
it's good to have a general feeling of the system by maximal localizing
the valence + conduction manifolds of MoS2.

Since we are focusing on the isolated manifolds of MoS2,
we only need to run [`max_localize`](@ref) function,
=#
U = max_localize(model);

# The initial spread is
omega(model)

# The final spread is
omega(model, U)

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
Random.seed!(12345)
R = Wannier.rand_unitary(eltype(model_val.U), size(model_val.U)...)
# and assign it to the `model_val`
model_val.U .= R;
# of course now we have totally destroyed the gauge 😈
omega(model_val)
# and we will ask `parallel_transport` to explicitly use our random matrices,
# by setting argument `use_U = true`,
U2, _ = parallel_transport(model_val; use_U=true)
#=
note the function returns a lit bit extra information which is irrelevant to this tutorial,
so we discard it by a `_`.

Then we inspect the new gauge
=#
omega(model_val, U2)

#=
## Fix gauge at first kpoint

However, the [`parallel_transport`](@ref) itself does not fix the gauge of the first kpoint,
which is arbitrary. We can fix it by maximal localizing w.r.t. a single (i.e. independent of kpoints)
`n_wann * n_wann` rotation matrix `W`, by calling [`opt_rotate`](@ref),
=#
model_val.U .= U2;
W = opt_rotate(model_val)
#=
This is convenient since the calculation is cheap, and helps evade local minimum.

Then let's rotate the input gauge by `W`,
=#
U3 = rotate_U(U2, W);
#=
and inspect the new gauge,
=#
omega(model_val, U3)

#=
## Final maximal localization

Often it is helpful to run a final maximal localization that could further
smoothen the gauge,
=#
model_val.U .= U3;
U4 = max_localize(model_val);
# and inspect the final gauge
omega(model_val, U4)

#=
## Band interpolation

Finally, let's have a look at the band interpolation.

Valence + conduction bands,
=#
model.U .= U;
interp_model = Wannier.InterpModel(model)
kpi, E = interpolate(interp_model)
# and top valence band,
model_val.U .= U4;
interp_model_val = Wannier.InterpModel(model_val)
kpi2, E2 = interpolate(interp_model_val)
# and plot the band structure
P = plot_band_diff(kpi, E, E2)
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
since we won't use the gauge from the previous `model.U`. We use
this `[1]` to specify that we only need 1 WF.

The initial spread is
=#
omega(model_top)
# with PTG,
U_top, _ = parallel_transport(model_top)
omega(model_top, U_top)
#=
In the single-band case, there is no need to run optimal rotation.
But we can still run maximal localization,
=#
model_top.U .= U_top;
U_top2 = max_localize(model_top)
omega(model_top, U_top2)
# and band interpolation
interp_model_top = Wannier.InterpModel(model_top)
kpi3, E3 = interpolate(interp_model_top)
# and compare
P = plot_band_diff(kpi, E, E3)
Main.HTMLPlot(P, 500)  # hide

#=
Again, rerun with a finer kgrid and see the improvements of interpolation quality.
Also, you can try to plot the WFs, by first truncating the `unk`
=#
Wannier.truncate_unk(".", [7], "truncate")
# the new `unk`s are stored in a subfolder `truncate`,
# now write the realspace WF
write_realspace_wf("wjl", model_top; n_supercells=[3, 3, 2], unkdir="truncate")
# and visualize with your favorite tool.

# As a reference, here is the real space WFs visualized with [`WannierPlots.jl`].
using JSServe  # hide
Page(; exportable=true, offline=true)  # hide
# First, load the plotting packages
using WGLMakie
set_theme!(; resolution=(800, 800))
using WannierPlots
#=
!!! tip

    Here we want to show the WFs in this web page, so we first load `WGLMakie`.
    When you use the `WannierPlots` package in REPL, you can first load `GLMakie`,
    then the WFs will be shown in a standalone window.
=#
# Read the 1st WF
xsf = read_xsf("wjl_00001.xsf");
# Visualize with `WannierPlots.jl`,
pos = inv(xsf.primvec) * xsf.atom_positions  # to fractional coordinates
atom_numbers = parse.(Int, xsf.atoms)  # to integer atomic numbers
plot_wf(xsf.rgrid, xsf.W, xsf.primvec, pos, atom_numbers)

#=
Now we have automated Wannierization of the valence manifold,
and the top valence band of MoS2, without any initial projection!
=#
