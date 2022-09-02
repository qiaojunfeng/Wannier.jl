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

## Fortran related

```@docs
isbinary
parse_float
parse_bool
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
rotate_gauge(O::Array{T,3}, A::Array{T,3}) where {T<:Number}
eyes_A
rotate_A
rotate_M
isunitary
get_projectability
findvector
```

## Kpoints

```@docs
get_kpoint_mappings
make_supercell
get_kpoints(kgrid::AbstractVector{Int}; fractional::Bool=true)
sort_kpoints
get_kgrid
```

## Centers

```@docs
find_nearests
find_nearest_atom
wrap_centers
```
