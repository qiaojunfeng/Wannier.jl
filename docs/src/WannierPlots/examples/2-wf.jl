# # 2. Realspace WFs of graphene

#=
```@meta
CurrentModule = WannierPlots
```
=#

#=
In this tutorial, we will plot realspace representation.

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`2-wf.ipynb`](./2-wf.ipynb)
=#

# ## Preparation
# Load the package
using Wannier: get_atom_number, read_xsf

using WGLMakie
set_theme!(; resolution=(800, 800))
using WannierPlots

#=
!!! tip

    Here we want to show the WFs in this web page, so we first load `WGLMakie`.
    When you use the `WannierPlots` package in REPL, you can first load `GLMakie`,
    then the WFs will be shown in a standalone window.
=#

using JSServe  # hide
Page(; exportable=true, offline=true)  # hide

# Path of current tutorial
PWD = "4-realspace"

# Read the 1st WF
xsf = read_xsf("$PWD/graphene_00001.xsf");
# Visualize with `WannierPlots.jl`,
pos = inv(xsf.primvec) * xsf.atom_positions  # to fractional coordinates
atom_numbers = get_atom_number(xsf.atoms)  # to integer atomic numbers
plot_wf(xsf.rgrid, xsf.W, xsf.primvec, pos, atom_numbers)
