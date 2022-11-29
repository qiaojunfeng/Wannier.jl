# Model

The [`Model`](@ref Model) `struct` abstracts a single Wannierization, it contains all the necessary
information for maximally localization of the crystal structure.

```@meta
CurrentModule = Wannier
```

## Contents

```@contents
Pages = ["model.md"]
Depth = 2
```

## Index

```@index
Pages = ["model.md"]
```

## Model struct

```@docs
Model
Model(lattice::Mat3{T}, atom_positions::Matrix{T}, atom_labels::Vector{String}, kgrid::Vec3{Int}, kpoints::Matrix{T}, bvectors::BVectors{T}, frozen_bands::AbstractMatrix{Bool}, M::Array{Complex{T}, 4}, U::Array{Complex{T}, 3}, E::Matrix{T}) where {T<:Real}
```

## Model functions

```@docs
rotate_gauge(model::Model, U::Array{T,3}) where {T<:Number}
truncate(model::Model, keep_bands::AbstractVector{Int})
```

## Spread

```@docs
Spread
omega
omega_grad
omega_local
center
position_op
berry_connection
```
