# Interpolation

Wannier interpolation starts by constructing one persistent
[`InterpolationModel`](@ref):

```julia
interpolation_model = InterpolationModel(
    model;
    real_space = MinimumDistance(),
)

result = interpolate(interpolation_model, kpoints, BandEnergy())
result.band_energy
```

`interpolate` always returns a named result. Every result array places the
k-point index on its final axis; `band_energy` therefore has shape
`n_wannier × n_kpoints` and is measured in eV.

For a labeled `KPath`, [`band_structure`](@ref) also retains the linear
path coordinate and high-symmetry labels:

```julia
bands = band_structure(interpolation_model, kpath)
bands.path_coordinate
bands.symmetry_point_labels
bands.band_energy
```

## Real-space scheme

[`MinimumDistance`](@ref) is the default and selects the nearest periodic image
for each pair of Wannier orbitals. [`WignerSeitz`](@ref) instead uses one common
Wigner--Seitz selection. Scheme-specific degeneracies and replica translations
are absorbed during construction; evaluation uses one ordinary Fourier sum on
the resulting [`RealSpaceDomain`](@ref).

Every primitive operator in an `InterpolationModel` uses the same domain. Its
coefficients have the layout

```text
n_wannier × n_wannier × component_shape... × n_Rvectors
```

so scalar, vector, and rank-two operators retain their natural physical axes.
Only numerical kernels flatten these axes.

## General operators

Additional lattice-periodic operators are supplied in Bloch space:

```julia
spin = BlochOperator(
    spin_matrices;
    law = AxialVector(time_reversal = Odd()),
    hermitian = true,
)

interpolation_model = InterpolationModel(
    model;
    operators = (; spin),
    real_space = MinimumDistance(),
)
```

Full input matrices have layout
`n_bands × n_bands × component_shape... × n_kpoints`. A diagonal scalar
operator may instead use `n_bands × n_kpoints`. The operator law and Hermiticity
are explicit so later symmetry closure can apply the correct transformation to
every matrix and physical-component index.

```@meta
CurrentModule = Wannier
```

## Public interface

```@autodocs
Modules = [Wannier]
Pages = [
    "interpolation/types.jl",
    "interpolation/observables/band_energy.jl",
    "interpolation/workflows/band_structure.jl",
]
```

## Migration status

The operator-specific `TBOperator` and interpolator types remain available while
the other observables are migrated to the common model. New code should use
`InterpolationModel` and `interpolate`.
