# Utility

These are some simple convenience functions for development.

```@meta
CurrentModule = Wannier
```

## Contents

```@contents
Pages = ["util.md"]
Depth = 2
```

## Index

```@index
Pages = ["util.md"]
```

## Atomic properties

```@docs
get_atom_number
```

## Lattice

```@docs
get_recip_lattice
get_lattice
```

## Matrices

```@docs
imaglog
orthonorm_lowdin
fix_global_phase
compute_imre_ratio
rotate_gauge(O::Array{T,3}, U::Array{T,3}) where {T<:Number}
eyes_U
rotate_U
rotate_M
isunitary
get_projectability
findvector
rand_unitary
```

## Kpoints

```@docs
get_kpoint_mappings
make_supercell
get_kpoints(kgrid::AbstractVector{Int}; fractional::Bool=true)
sort_kpoints
get_kgrid
```

## KPath

The `KPath` and `KPathInterpolant` are defined in `Brillouin.jl`,
they are used to store the high-symmetry kpoints and their labels.

```@docs
KPath(
    lattice::AbstractMatrix, kpoint_path::Vector{Vector{Pair{Symbol,T}}}
) where {T<:AbstractVector{<:Real}}
interpolate_w90
get_x
get_kpoints(kpi::KPathInterpolant)
get_kpath
```

## Centers

```@docs
find_nearests
find_nearest_atom
wrap_centers
```
