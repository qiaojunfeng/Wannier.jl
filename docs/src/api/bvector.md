# K-space stencil

The b vectors connect each source k point to neighboring points used to evaluate
Wannier centers, spreads, and their gradients. [`KSpaceStencil`](@ref) stores
this finite-difference geometry over the complete source k-space grid; it is not
an evaluation grid for interpolation.

[`generate_kspace_stencil`](@ref) groups candidate b vectors into equal-distance
shells, removes redundant parallel shells, solves the completeness condition,
and reproduces Wannier90's ordering at every k point. The `atol` keyword controls
the floating-point comparisons and corresponds to Wannier90's `kmesh_tol`.

```@meta
CurrentModule = Wannier
```

```@autodocs
Modules = [Wannier]
Pages = [
    "kpoints/kstencil_shell.jl",
    "kpoints/kstencil.jl",
]
```
