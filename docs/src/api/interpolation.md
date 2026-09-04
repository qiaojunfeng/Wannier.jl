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

Several observables can be requested together:

```julia
result = interpolate(
    interpolation_model,
    kpoints,
    (BandEnergy(), BandVelocity()),
)
result.band_energy
result.band_velocity
```

The combined calculation shares its Fourier phase block, Hamiltonian
derivatives, and eigendecomposition. `band_velocity` stores
``\hbar \mathbf{v}_n = \partial\varepsilon_n/\partial\mathbf{k}`` in eV Å,
with shape `3 × n_wannier × n_kpoints`; the leading axis contains Cartesian
components.

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

## Symmetry closure

Construct [`WannierSymmetry`](@ref) from standardized `.isym` data and the
prescribed fractional Wannier centers, then pass it to the interpolation-model
constructor:

```julia
symmetry = WannierSymmetry(isym, fractional_centers)
interpolation_model = InterpolationModel(
    model;
    real_space = MinimumDistance(),
    symmetry,
)
```

Construction completes the real-space support under unitary, antiunitary, and
Hermitian partners and applies the symmetry projector once. Subsequent queries
remain ordinary Fourier sums while satisfying off-mesh covariance. Scalar,
polar-vector, axial-vector, and general Cartesian-tensor laws are supported;
their time-reversal parity is part of the operator law. Minimum-distance ties
are transported with their orbital-pair-dependent cell shifts and weights.

`SymmetryConstraint` used by SAWF localization contains this same
`WannierSymmetry` as `constraint.wannier_symmetry`; its remaining tables are
specific to the IBZ optimization.

```@docs
Wannier.WannierSymmetry
```

## Wannier90 real-space data

Already constructed Wannier90 Hamiltonians enter the same model without a
legacy tight-binding wrapper. Supply the parsed `HrDat` or `TbDat`, fractional
Wannier centers, and the matching `WsvecDat` when it is available:

```julia
hrdat = read_w90_hr_dat(prefix * "_hr.dat")
wsvec = read_w90_wsvec_dat(prefix * "_wsvec.dat")
interpolation_model = InterpolationModel(
    hrdat,
    lattice;
    fractional_centers,
    wsvec,
)
```

The file degeneracies and minimum-distance translations are absorbed directly
into the common `RealSpaceDomain` representation during construction.

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

Wannier90 spin matrices have a direct adapter:

```julia
spin = BlochOperator(read_spn(prefix * ".spn"))
```

Full input matrices have layout
`n_bands × n_bands × component_shape... × n_kpoints`. A diagonal scalar
operator may instead use `n_bands × n_kpoints`. The operator law and Hermiticity
are explicit so symmetry closure can apply the correct transformation to
every matrix and physical-component index.

The spin primitive is consumed by [`SpinExpectation`](@ref). Without an axis,
the result has shape `3 × n_wannier × n_kpoints`; supplying a Cartesian axis
returns its projection with shape `n_wannier × n_kpoints`:

```julia
result = interpolate(
    interpolation_model,
    kpoints,
    (
        BandEnergy(),
        SpinExpectation(),
    ),
)

spin_z = interpolate(
    interpolation_model,
    kpoints,
    SpinExpectation([0, 0, 1]),
)
```

Both requests reuse the Hamiltonian eigensystem. Spin values inherit the units
of the primitive operator and are dimensionless for Wannier90 `spn` data.

The finite-difference Wannier-gauge Berry connection is constructed directly
from the model overlaps and gauges:

```julia
interpolation_model = InterpolationModel(
    model;
    operators = (;
        berry_connection = BerryConnection(; imlog_diag = false),
    ),
    real_space = MinimumDistance(),
)
```

`BerryConnection` is a construction recipe rather than a `BlochOperator`: a
connection has an inhomogeneous gauge-transformation law and cannot be treated
as a lattice-periodic Cartesian vector. Its coefficients use the physical
layout `n_wannier × n_wannier × 3 × n_Rvectors`. When symmetry is supplied,
construction subtracts the prescribed Wannier-center term, closes and projects
the resulting homogeneous polar vector, and restores the affine center term.
The interpolated connection therefore obeys its proper affine law rather than
being incorrectly projected as an ordinary vector.

Berry curvature is selected as an observable formulation. The default
occupied-manifold recipe uses the Lopez--Vanderbilt--Thonhauser--Souza formula
and returns an antisymmetric Cartesian tensor of shape `3 × 3 × n_kpoints` in
Å²:

```julia
occupied = interpolate(
    interpolation_model,
    kpoints,
    BerryCurvature(fermi_energy),
)

by_band = interpolate(
    interpolation_model,
    kpoints,
    BerryCurvature(; formulation = WYSV06BandResolved()),
)
```

The band-resolved result has shape
`3 × 3 × n_wannier × n_kpoints`. `WYSV06()` selects the alternative
occupied-manifold formulation. Berry curvature shares the Hamiltonian,
Fourier phases, eigensystem, and the requested Hamiltonian and connection
derivatives with other observables in the same call.

Orbital magnetization additionally needs the two Hamiltonian-weighted position
moments derived from a Wannier90 `uHu` file:

```julia
uHu = read_uHu(prefix * ".uHu").uHu
interpolation_model = InterpolationModel(
    model;
    operators = (;
        berry_connection = BerryConnection(; imlog_diag = false),
        hamiltonian_position = HamiltonianPosition(),
        position_hamiltonian_position = PositionHamiltonianPosition(uHu),
    ),
    real_space = MinimumDistance(),
)

result = interpolate(
    interpolation_model,
    kpoints,
    OrbitalMagnetization(fermi_energy),
)
```

`result.orbital_magnetization` is the antisymmetric Cartesian integrand with
shape `3 × 3 × n_kpoints` and units eV Å². Under symmetry, the two moment
recipes are centered on the Wannier functions adjacent to their coordinate
legs. Their centered polar-vector and rank-two-tensor coefficients are closed
and projected jointly with the Hamiltonian before the conventional uncentered
moments are restored.

```@meta
CurrentModule = Wannier
```

## Public interface

```@autodocs
Modules = [Wannier]
Pages = [
    "interpolation/types.jl",
    "io/w90/interpolation.jl",
    "interpolation/planning.jl",
    "interpolation/observables/band_energy.jl",
    "interpolation/observables/band_velocity.jl",
    "interpolation/observables/spin_expectation.jl",
    "interpolation/operators/berry_connection.jl",
    "interpolation/operators/orbital_magnetization.jl",
    "interpolation/observables/berry_curvature.jl",
    "interpolation/observables/orbital_magnetization.jl",
    "interpolation/workflows/band_structure.jl",
]
```

## Migration status

The operator-specific `TBOperator` and interpolator types remain temporarily as
reference paths while callers are migrated. All currently supported pointwise
physical quantities are available through `InterpolationModel` and
`interpolate`; new code should use this interface.
