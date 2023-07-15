# # Introduction to Julia artifacts

#=
```@meta
CurrentModule = Wannier
```
=#

#=
Throughout the examples, we need to load input files (e.g., wannier90
`amn`/`mmn`/`eig` files) to construct `Model`s for Wannierizations.
Often, preparing these input files requires density-functional-theory
calculations which can be time-consuming, and distracting ourselves from
the main purpose of Wannierizations and Wannier interpolations.

Julia has a convenient [`artifact`](https://pkgdocs.julialang.org/v1/artifacts/)
system allowing us to load binaries or data files easily in Julia script or REPL,
without the need of manually downloading and placing them in the right folder.

In `Wannier.jl`, we use such system to load our pre-computed
[`WannierDatasets`](https://github.com/qiaojunfeng/WannierDatasets), which
contains some typical materials, e.g., silicon, copper, graphene, etc.
Thus, we can easily load these datasets and construct WFs hassle-free.
=#

# ## How to use `artifact`
using LazyArtifacts

using FileTrees
FileTree(artifact"Si2")
