# # 6. Splitting the valence/conduction of silicon

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In the previous tutorial, we have seen how to use the `parallel_transport` function
to automate the construction of WFs for isolated manifold.
In this tutorial, we will again use this technique, combined with disentanglement,
to Wannierize the valence and conduction manifolds of silicon, separately.

Usually, Wannierization of valence or conduction manifold is difficult because
there is no prior knowledge of the bonding or anti-bonding orbitals, except for
simple system.

For the valence manifold, since it is isolated, we can

- use random Gaussian as initial projections, then run many iterations of maxiaml localization,
- use SCDM for the isolated manifold, then run maximal localization.

However, for the conduction manifold, it is entangled with higher energy bands,
so the random Gaussian is not applicable. In such case, we can

- use SCDM with specified ``\mu`` and ``\sigma`` for the conduction manifold,
    then run maximal localization.

However, the quality of Wannierization sensitively depends on the ``\mu`` and ``\sigma``,
and also SCDM works in real space which is more memory consuming.

Now, our solution is

- run the normal disentanglement on valence + conduction manifold, the initial guess
    can be the analytic spdf projectors, or using numerical atomic orbitals (NAOs)
    from pseudopotentials (for plane-wave code like QE, this automates the initial
    projections for the disentanglement).
- separate the valence + conduction manifold into two, by diagonalizing the Wannier
    gauge Hamiltonian, the eigenvectors are unitary matrices, in which the upper (lower)
    part maps to the valence (conduction) manifold. Thus we have two separated submanifolds.
- run the `parallel_transport` on the separated submanifolds, to construct the
    parallel transport gauge (PTG) for valence and conduction, respectively.
    We can further run optimal rotation and maximal localization to smoothen the gauge.

The disentanglement generates a smooth "total" manifold, helping us getting rid of
higher energy states, and the PTG automates the gauge construction for the two separated
submanifolds. In final, we have an automated approach to Wannierize the valence and/or
conduction manifolds without human intervention.

## Outline

1. construct a [`Model`](@ref), by reading the `win`, `amn`, `mmn`, and `eig` files
2. split the `Model` by diagonalizing the Wannier gauge Hamiltonian
3. run [`parallel_transport`](@ref) on the two `Model`s to construct PTG for valence
    and conduction, respectively
4. maxiaml localize to smoothen the PTG
5. interpolate band structures

!!! note

    These tutorials assume you have already been familiar with the
    Wannierization workflow using `QE` and `Wannier90`, a good starting
    point could be the tutorials of
    [`Wannier90`](https://github.com/wannier-developers/wannier90).

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`6-split.ipynb`](./6-split.ipynb)
    - Julia script: [`6-split.jl`](./6-split.jl)
=#

# ## Preparation
# Load the package
using Wannier
using WannierPlots

# Path of current tutorial
CUR_DIR = "6-split"

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = read_w90("$CUR_DIR/si2")

#=
## Disentanglement

First let's disentangle the valence + conduction manifold,
=#
A = disentangle(model);

#=
## Splitting the model

Then we split the model into two, by diagonalizing the Wannier gauge Hamiltonian,
we need to pass it the number of valence bands, `n_val`,
=#
n_val = 4
model_v, model_c, Uv, Uc = split_model(model, n_val)
#=
Returned are two separated `Model`s, and the two corresponding gauge transformation
matrices from the total manifold to the two separated submanifolds.
This is useful if you want to further rotate the gauge of some operators,
or if you want to rotate the `UNK` files for plotting WFs,
in that case you need to pass these two matrices to the function
[`split_unk`](@ref).

Next, we construct PTG for the two `Model`s,
=#
Av, _ = parallel_transport(model_v)
model_v.A .= Av;
# and
Ac, _ = parallel_transport(model_c)
model_c.A .= Ac;
# and take a look at the spread
omega(model_v)
omega(model_c)
#=
!!! tip

    There is a [`split_wannierize`](@ref) function that combines the above steps,
    however, for pedagogical purpose we run them manually.

## Further maximal localization

We can further run maximal localization to smoothen the gauge,
=#
Av = max_localize(model_v)
model_v.A .= Av;
# and
Ac = max_localize(model_c)
model_c.A .= Ac;
#=
!!! tip

    For such simple case (few bands, very little kpoints) a direct maximal localization
    is sufficient. But for more complex scenarios, running an [`opt_rotate`](@ref) might
    be helpful.

## Comparing spreads

The valence + conduction,
=#
model.A .= A;
omega(model)
# the valence,
omega(model_v)
# the conduction,
omega(model_c)
#=
Look, we have well-localized WFs for all the three cases,
the initial 8 s,p orbitals are separated into 4 bonding
and 4 anti-bonding orbitals! 🚀

## Band structure

Finally, let's compare band interpolation,
=#
using PlotlyJS
# the valence + conduction,
interp_model = Wannier.InterpolationModel(model)
kpi, E = interpolate(interp_model)
# the valence,
interp_model_v = Wannier.InterpolationModel(model_v)
_, Ev = interpolate(interp_model_v)
# the conduction,
interp_model_c = Wannier.InterpolationModel(model_c)
_, Ec = interpolate(interp_model_c)
# and compare valence,
P = plot_band_diff(kpi, E, Ev)
Main.HTMLPlot(P, 500)  # hide
# and conduction,
P = plot_band_diff(kpi, E, Ec)
Main.HTMLPlot(P, 500)  # hide
#=
!!! note

    Again, the interpolation quality is not good (especially for the conduction
    which has much larger spread WFs), bacause of the very coarse
    kgrid, increase the kgrid and compare the resulting band structures.

Bravo, we have successfully Wannierized the valence and conduction manifolds,
by splitting the gauge of a valence + conduction calculation!
=#
