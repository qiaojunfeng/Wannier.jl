# Wannierize

These are some Wannierization algorithms operating on [`Model`](@ref Model).

```@meta
CurrentModule = Wannier
```

## Contents

```@contents
Pages = ["wannierize.md"]
Depth = 2
```

## Index

```@index
Pages = ["wannierize.md"]
```

## Disentanglement

!!! note

    In general, the user only need to call `set_frozen_win!` to set frozen window for the `Model`,
    then call `disentangle` to disentangle the `Model`. Other functions are for internal use.

```@docs
get_frozen_bands
set_frozen_degen!
check_frozen_bands
set_frozen_win!
get_frozen_proj
set_frozen_proj!
orthonorm_freeze
X_Y_to_U
U_to_X_Y
XY_to_X_Y
X_Y_to_XY
omega(bvectors::BVectors{FT},M::Array{Complex{FT},4},X::Array{Complex{FT},3},Y::Array{Complex{FT},3}) where {FT<:Real}
omega_grad(bvectors::BVectors{FT},M::Array{Complex{FT},4},X::Array{Complex{FT},3},Y::Array{Complex{FT},3},frozen::BitMatrix) where {FT<:Real}
zero_froz_grad!
get_fg!_disentangle
disentangle
```

## Maximal localization

!!! note

    In general, the user only need to call `max_localize` to maximally localize the `Model`.
    Other functions are for internal use.

```@docs
get_fg!_maxloc
max_localize
```

## Parallel transport

!!! note

    In general, the user only need to call `parallel_transport` to construct the parallel transport gauge
    for the `Model`.

```@docs
parallel_transport
compute_error
```

## Optimal rotation

!!! note

    In general, the user only need to call `opt_rotate` to find the optimal rotation matrix
    for the `Model`, then call `rotate_U` to rotate the initial `U` matrices at each kpoint.

```@docs
get_fg!_rotate
opt_rotate
rotate_U(U::Array{T,3}, W::Matrix{T}) where {T<:Complex}
```

## Splitting the Model

!!! note

    In general, the user only need to call `split_wannierize` to construct two `Model`s
    for valence and conduction, respectively.

```@docs
split_eig
split_unk
split_model
split_wannierize
```

## Constraining WF centers

### Disentanglement

!!! note

    In general, the user only need to call `disentangle_center` to disentangle with penalty
    on WF centers.

```@docs
omega_center(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
omega_center_grad(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
    frozen::BitMatrix,
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
get_fg!_center_disentangle
disentangle_center
```

### Maximal localization

!!! note

    In general, the user only need to call `max_localize_center` to disentangle with penalty
    on WF centers.

```@docs
SpreadCenter
omega_center(
    bvectors::BVectors{T},
    M::Array{Complex{T},4},
    U::Array{Complex{T},3},
    r₀::Matrix{T},
    λ::T,
) where {T<:Real}
omega_center_grad(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    U::Array{Complex{FT},3},
    r::Matrix{FT},
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
omega_center_grad(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    U::Array{Complex{FT},3},
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
get_fg!_center_maxloc
max_localize_center
```
