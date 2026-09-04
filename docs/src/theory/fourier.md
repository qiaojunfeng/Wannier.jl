# Wannier interpolation

```@meta
CurrentModule = Wannier
```

Wannier interpolation represents a smooth, lattice-periodic Bloch-space
operator by a finite set of real-space matrix elements. In Wannier.jl,
[`InterpolationModel`](@ref) owns these matrix elements and [`interpolate`](@ref)
is the single pointwise evaluation interface.

## From the sampled mesh to real space

Let ``O^{\mathrm B}(\mathbf{k})`` be an operator on the uniform source mesh and
``U(\mathbf{k})`` the Wannier gauge. Its Wannier-gauge matrix is

```math
O^{\mathrm W}(\mathbf{k}) =
U^\dagger(\mathbf{k}) O^{\mathrm B}(\mathbf{k}) U(\mathbf{k}).
```

For an operator diagonal in the input eigenstate basis, such as the
Hamiltonian, the middle factor is simply the diagonal matrix of sampled
eigenvalues. The quotient-lattice coefficients are

```math
\widehat O^{\mathrm W}_{mn}(\mathbf R) =
\frac{1}{N_k}\sum_{\mathbf k}
e^{-2\pi i\,\mathbf k\cdot\mathbf R}
O^{\mathrm W}_{mn}(\mathbf k).
```

Here both ``\mathbf k`` and the integer lattice vector ``\mathbf R`` use
fractional coordinates. Vectors that differ by a vector of the real-space
superlattice dual to the sampling mesh represent the same discrete Fourier
coefficient, so a real-space selection rule chooses finite representatives and
assigns their weights.

[`WignerSeitz`](@ref) selects representatives in the Wigner--Seitz cell of the
sampling superlattice and divides boundary coefficients among equidistant
representatives. [`MinimumDistance`](@ref), the default, selects the nearest
periodic representatives separately for every pair of Wannier centers. The
latter usually reduces phase errors for off-mesh interpolation. Both schemes
are construction choices:

```julia
ws_model = InterpolationModel(model; real_space = WignerSeitz())
md_model = InterpolationModel(model; real_space = MinimumDistance())
```

All primitive operators in an interpolation model are remapped onto one common
finite [`Wannier.RealSpaceDomain`](@ref). Callers therefore do not manipulate
scheme-specific real-space containers or degeneracy tables.

## Evaluation

After representative selection, every primitive operator has coefficients
``O^{\mathrm W}_{mn}(\mathbf R)`` and is evaluated by the ordinary Fourier sum

```math
O^{\mathrm W}_{mn}(\mathbf k) =
\sum_{\mathbf R} e^{2\pi i\,\mathbf k\cdot\mathbf R}
O^{\mathrm W}_{mn}(\mathbf R).
```

One phase matrix is shared by all operators requested in the same call. The
Hamiltonian is diagonalized once, and observables reuse its eigensystem. For
example,

```julia
result = interpolate(
    interpolation_model,
    kpoints,
    (BandEnergy(), BandVelocity(), SpinExpectation()),
)
```

Cartesian derivatives do not require a separate interpolation object. They
follow by differentiating the phase,

```math
\frac{\partial O^{\mathrm W}(\mathbf k)}{\partial k_\alpha}
= i\sum_{\mathbf R} R_\alpha
e^{i\mathbf k\cdot\mathbf R} O^{\mathrm W}(\mathbf R),
```

where the last expression uses Cartesian ``\mathbf k`` and ``\mathbf R``.
The planner computes these derivative sums only for observables that need them.

## Symmetry-preserving construction

Supplying `symmetry` to [`InterpolationModel`](@ref) is a strong guarantee. At
construction time, the real-space support is closed under the symmetry group
and Hermitian conjugation, and each coefficient orbit is projected according
to the operator's scalar, vector, or tensor transformation law. Ordinary
Fourier evaluation then obeys the requested covariance at arbitrary off-mesh
points; no group average is performed for each query.

The Hamiltonian is only the time-reversal-even scalar instance of this general
operator construction. Position, spin, and Hamiltonian-weighted moments use the
same real-space domain and evaluation kernel, while their observable recipes
assemble Berry curvature, orbital magnetization, and other derived quantities.
