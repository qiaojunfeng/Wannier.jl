# # 4. Realspace WFs of graphene

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will disentangle the 2D graphene, and compute its realspace representation.

1. construct a [`Model`](@ref) for `Wannier.jl`, by reading the `win`, `amn`, `mmn`, and `eig` files
2. run `Wannier.jl` [`disentangle`](@ref) on the `Model` to minimize the spread
3. write the real space WFs to `xsf` files

!!! note

    These tutorials assume you have already been familiar with the
    Wannierization workflow using `QE` and `Wannier90`, a good starting
    point could be the tutorials of
    [`Wannier90`](https://github.com/wannier-developers/wannier90).

    Different from previous tutorials, the `amn/mmn/eig/unk` files in this tutorial
    are generated by the `DFTK.jl` package. See `dftk.jl` script in the tutorial folder
    for more information.

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`4-realspace.ipynb`](./4-realspace.ipynb)
    - Julia script: [`4-realspace.jl`](./4-realspace.jl)
=#

# ## Preparation
# Load the package
using Wannier

# Path of current tutorial
CUR_DIR = "4-realspace"

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = read_w90("$CUR_DIR/graphene")

#=
## Disentanglement and maximal localization

The [`disentangle`](@ref) function
will disentangle and maximally localize the spread
functional, and returns the gauge matrices `A`,
=#
A = disentangle(model)

# The initial spread is
omega(model)

# The final spread is
omega(model, A)

#=
!!! note

    The convergence thresholds is determined by the
    keyword arguments of [`disentangle`](@ref), e.g., `f_tol` for the tolerance on spread,
    and `g_tol` for the tolerance on the norm of spread gradient, etc. You can use stricter thresholds
    to further minimize a bit the spread.
=#

#=
## Write real space WFs

Now save the `A` to the `model`,
=#
model.A .= A;
#=
!!! tip

    You need to use `.=` to assign the `A` to the `model`,
    because in Julia `model` is an immutable `struct` so you cannot
    use `model.A = A`.
=#

#=
The [`write_realspace_wf`](@ref) function reads the `UNK` files,
compute the real space WFs in a `n_supercells`-sized super cell,
and write them to `xsf` files,
=#
write_realspace_wf("$CUR_DIR/wjl", model; unkdir=CUR_DIR, n_supercells=3, format=:xsf)

#=
Now, open the `wjl_00001.xsf`, etc. files with a 3D
visualizer, e.g., `vesta`, to have a look at the WFs!

!!! note

    To save space, we use a very coarse kpoint grid and
    reduced sampling of the real space grid, thus the resolution
    are bad, rerun the calculation with a finer kpoint grid and
    see the beautiful WFs 😎

=#

#=
## Compute WF centers in realspace

There are some other functions that might be useful for evaluating operators
in real space, e.g., computing WF centers.

First we need to read the `UNK` files, and construct the real space WFs
in a `3 * 3 * 1`-sized super cell (i.e., `model.kgrid`),
=#
rgrid, W = read_realspace_wf(model, model.kgrid, CUR_DIR)
#=
The real space WFs `W`, are defined on the grid `rgrid`.

To compute WF centers, invoke
=#
center(rgrid, W)
# columns are the WF centers in Cartesian coordinates.

# Compare with WF center computed in reciprocal space,
center(model)

#=
and yes, they are different, the z coordinate is wrongly computed
in the real space because the WFs are truncated along z (you can
see this by a 3D visualizer). If we translate the `rgrid` by
half of the `c` axis along z, then the WFs are complete, and the
real space WF centers are much closer to the reciprocal space results.

!!! note

    Because we are using a coarsely sampled real space grid in
    the `UNK` files, the WF centers can be improved by
    rerun the calculation with a more refined real space grid.
=#

#=
That's all about real space!
=#