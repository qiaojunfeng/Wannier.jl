# B vector

The ``\bm{b}``-vectors connect kpoint ``\bm{k}`` to its neighboring kpoints ``\bm{k}+\bm{b}``,
for calculating WF centers, spreads in reciprocal space.

The bvectors are first arranged in layers of shells, which contain bvectors having same norm.
Then parallel shells are deleted, and the shells satisfying B1 condition are the final bvectors.
At last, the bvectors are sorted at each kpoint, to reproduce exactly the same order as `Wannier90`.
This ensures that we have the same order as `mmn` file.

!!! note

    To reproduce the same order as `Wannier90`, we need to be careful with some floating point
    comparison, i.e., the `atol` keyword arguments in the following functions. Their default
    value reproduce the `Wannier90` order.
    Note the `kmesh_tol` in the `win` file also influence the search of bvectors.

!!! tip
    In most cases, the user only need to call [`get_bvectors`](@ref get_bvectors) to generate
    bvectors having the same order as `Wannier90`. Other functions are intermediate steps that are
    called inside `get_bvectors`.

```@meta
CurrentModule = Wannier
```

## Contents

```@contents
Pages = ["bvector.md"]
Depth = 2
```

## Index

```@index
Pages = ["bvector.md"]
```

## B vector shells and B vectors

```@autodocs
Modules = [Wannier]
Pages   = ["bvector.jl"]
```
