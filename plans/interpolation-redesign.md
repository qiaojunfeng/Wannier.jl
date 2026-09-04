# Wannier interpolation redesign

## Status and decisions

This plan replaces the current interpolation interface without preserving backward
compatibility. The redesign is allowed to remove or rename exported types and
functions.

The following interface decisions are fixed:

- `interpolate` is the single pointwise interpolation verb.
- `band_structure` remains a higher-level workflow for a path with labels and path
  coordinates.
- Observable recipe types use singular names: `BandEnergy`, `BerryCurvature`,
  `BandVelocity`, `SpinExpectation`, and `OrbitalMagnetization`.
- Result fields also use singular names: `result.band_energy`,
  `result.berry_curvature`, and so on.
- The public keyword is `real_space`, with scheme values `WignerSeitz()` and
  `MinimumDistance()`. "Real-space realization" can remain mathematical
  terminology in the documentation, but `Realization` is not part of the public
  interface.
- The persistent data type is `InterpolationModel`; there is no separate
  `Interpolator` wrapper.
- `KspaceStencil` and `KspaceStencilShells` are renamed to `KSpaceStencil` and
  `KSpaceStencilShells`.
- One `RealSpaceDomain` is shared by every primitive operator in an
  `InterpolationModel`.
- Supplying symmetry means that the constructed real-space operators have
  symmetry-closed support and symmetry-projected coefficients. There is no public
  Boolean that merely promises symmetry.
- Physical component axes retain their logical tensor shape in stored arrays and
  results. They are flattened only inside numerical kernels.
- A new physical quantity does not add a field or a subtype of
  `InterpolationModel`.
  It adds an observable recipe and, only when physically necessary, a primitive
  real-space operator.

## Implementation progress

As of September 2026, Phase 1 is complete. `InterpolationModel` constructs the
Hamiltonian through the general `BlochOperator` path, both real-space schemes
produce packed operators on one common `RealSpaceDomain`, and ordinary
`BandEnergy`/`band_structure` evaluation is blocked over k points. Wannier90
`HrDat` and `TbDat` file data enter through direct construction Adapters rather
than legacy R-space or interpolator objects. The legacy interpolation types
remain only for observables awaiting migration in Phases 3--5.

Phase 2 is also complete for the migrated observables. A typed internal plan
unions observable requirements, reuses one phase block and one Hermitian
eigendecomposition sweep, and requests Cartesian Hamiltonian derivatives only
when needed. Its allocation-reusing `interpolate!` path has fixed per-batch
bookkeeping (measured at 1200 bytes for batches from 1 to 64 k points), rather
than per-k-point temporary arrays. On the `Si2_valence` reference mesh, a
combined `BandEnergy`/`BandVelocity` request is about 1.8 times faster than two
separate requests in a single-threaded warm benchmark.

Phase 3's reusable symmetry and closure machinery is implemented.
`WannierSymmetry` now owns only basis-intrinsic operations, centers, orbital
representations, and cell shifts, and `SymmetryConstraint` references it while
retaining localization-specific tables. Construction closes and projects
Hermitian scalar, polar-vector, axial-vector, and arbitrary-rank Cartesian
tensor operators, including antiunitary parity and minimum-distance replicas.
The ordinary Fourier evaluation of closed Si2 HSE operators agrees with an
independent per-query projector and satisfies random off-mesh covariance at
roughly `3e-14`; sampled-mesh changes remain at the cleaned `.isym` data floor.
The symmetry-closed Si2 Hamiltonian also recovers its two cubic triplets at
Gamma to `1e-14`, closing the Phase 3 exit criterion. The planned graphene
dataset remains a separate end-to-end interpolation demonstration.

Phase 4 now includes `SpinExpectation`, `BerryConnection`, `BerryCurvature`,
and `OrbitalMagnetization`. The planner evaluates arbitrary required primitive
operators and their requested first derivatives on its shared phase block,
while band energy, velocity, spin, Berry curvature, and orbital magnetization
reuse one Hamiltonian eigensystem.
A direct `BlochOperator(::WannierIO.Spn)` adapter supplies the axial,
time-reversal-odd spin primitive. The dedicated Berry-connection construction
recipe preserves its logical vector axis without pretending that a connection
has a homogeneous pointwise operator law. WYSV06, band-resolved WYSV06, and
LVTS12 curvature recipes match the existing Fe SOC/postw90 references.
Hamiltonian-weighted position-moment construction recipes and the orbital
magnetization integrand also match that reference. Connection and
Hamiltonian-weighted moment closure now use the centered polar-vector and
rank-two-tensor laws derived in the manuscript; random off-mesh tests satisfy
covariance at `1e-11` or better. Migrating remaining callers and deleting the
legacy observable-specific interpolators remain for Phase 5.

Phase 0 numerical references are covered by the existing interpolation tests
and benchmark entry point; finer per-stage timing remains to be recorded. The
next implementation step is Phase 5 caller migration and deletion.

## Goals

1. Present one small, coherent Interface for Hamiltonian and general-operator
   interpolation.
2. Treat the Hamiltonian as a scalar, time-reversal-even instance of the general
   operator construction rather than as an independent Fourier implementation.
3. Support Wigner--Seitz and minimum-distance real-space schemes without exposing
   intermediate R-space storage types or a `simplify` operation.
4. Guarantee exact off-mesh covariance when symmetry is requested by closing and
   projecting real-space coefficient orbits once during construction.
5. Share Fourier phases, Hamiltonian derivatives, eigendecompositions, and gauge
   rotations across simultaneously requested observables.
6. Use packed, strided storage suitable for BLAS and future accelerator kernels.
7. Make extension with a new general operator or observable local: the core
   `InterpolationModel` definition and unrelated observables must not change.

## Non-goals

- Preserve the existing `TBOperator`, `BareRspace`, `MDRSRspace`,
  `HamiltonianInterpolator`, or other observable-specific interpolator interfaces.
- Expose Fourier-transform details, symmetry-orbit construction, dependency-graph
  nodes, or workspace management in the initial public Interface.
- Add GPU execution in the first implementation. The storage and kernel shapes
  must not preclude it.
- Change the physical formulas used for existing Berry-curvature and
  orbital-magnetization implementations except where tests reveal a defect.

## Target public Interface

### Ordinary band interpolation

```julia
interpolation_model = InterpolationModel(model; real_space = MinimumDistance())

result = interpolate(interpolation_model, kpoints, BandEnergy())
result.band_energy

bands = band_structure(interpolation_model, kpath)
```

`InterpolationModel(model, ...)` automatically constructs the Hamiltonian from
the model eigenvalues and gauge. Internally this follows the same `BlochOperator`
path as every other lattice-periodic operator.

### Symmetry-preserving interpolation

```julia
interpolation_model = InterpolationModel(
    model;
    real_space = MinimumDistance(),
    symmetry = wannier_symmetry,
)
```

The presence of `symmetry` is a strong invariant: every stored operator has a
declared transformation law, its support is closed under that law and Hermitian
conjugation, and its coefficients have been projected onto the corresponding
covariant subspace. Ordinary Fourier evaluation must therefore satisfy covariance
at arbitrary k points without a per-query group average.

### General operators

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
    symmetry = wannier_symmetry,
)
```

The `operators` keyword accepts a `NamedTuple`. Its keys become stable operator
identities and are also available for inspection. The values describe input
Bloch-space matrix elements and their mathematical transformation laws; users do
not construct real-space coefficient tables manually in the common path.

When `symmetry` is present, a user-supplied operator without a transformation law
is an error. The implementation must never silently assume scalar or time-even
behavior.

### One or several observables

```julia
result = interpolate(
    interpolation_model,
    kpoints,
    (
        BandEnergy(),
        BandVelocity(),
        SpinExpectation(),
        BerryCurvature(),
    ),
)

result.band_energy
result.band_velocity
result.spin_expectation
result.berry_curvature
```

`interpolate` always returns a named result, including for one observable. This
keeps the return convention stable as a calculation grows from one quantity to
several. Every result array uses the k-point index as its final axis. Each
observable documents all preceding axes and its physical units.

### Repeated high-performance evaluation

The first implementation keeps planning internal so the public Interface stays
deep. `interpolate` builds a typed execution plan from the observable tuple and
uses internal workspaces. After profiling real applications, an explicit
`plan_interpolation` Interface may be exposed if repeated plan construction is
measurable. Do not expose it preemptively.

An internal allocation-reusing form is required:

```julia
Wannier.interpolate!(result_view, plan, kpoint_batch, workspace)
```

It is qualified and unexported initially. It operates on one k-point batch and a
matching view of the destination; it must not allocate or retain arrays sized by
the total number of requested k points. The allocating public method delegates to
it.

## Dense-grid readiness contract

Million-point grids, streamed output, symmetry-reduced sampling, and FFT execution
are implemented in the follow-on plan
[`dense-grid-interpolation.md`](dense-grid-interpolation.md). The present redesign
does not expose those additional Interface elements, but it must establish the
following invariants so that the follow-on work replaces no core data structure or
observable formula.

1. **No full-mesh intermediates.** Only a final collected result may scale with
   the total number of requested k points. Fourier phases, primitive operator
   values, derivatives, eigensystems, and observable intermediates scale with an
   internal batch size.
2. **Abstract k-point input.** `interpolate` accepts an indexable
   `AbstractVector`-like collection and accesses it by ranges. It never calls
   `collect` on the complete input. A lazy regular-grid type can therefore satisfy
   the same Interface later.
3. **Batch-oriented implementation.** The internal flow is k-point batch to
   primitive-operator batch to shared eigensystem/intermediates to observable
   batch to destination view. Observable formulas do not own Fourier loops.
4. **Destination-independent assembly.** Observable assembly writes into result
   views for a supplied k-point range. It does not assume that the complete result
   is an in-memory `Array`.
5. **Final k-point axis.** Every result places the k-point index last, so one batch
   maps to a natural slice such as `result.berry_curvature[:, :, k_range]`.
6. **Arbitrary real-space domain.** Fourier kernels make no rectangular-grid or
   fundamental-cell assumption about `RealSpaceDomain`; minimum-distance and
   symmetry-generated R vectors remain valid.
7. **Separated primitive evaluation and assembly.** The direct Fourier kernel
   produces a primitive-operator batch consumed by observable assembly. A later
   FFT implementation can produce the same batch representation without changing
   any observable recipe.
8. **Configurable batching.** Batch size is an internal plan parameter rather than
   a constant or a function of the total k-point count. The follow-on plan will
   derive it from an explicit memory limit.

These are implementation invariants, not additional public concepts. They deepen
the interpolation Module by preserving one observable Interface across small
paths, arbitrary point clouds, and future dense-grid execution.

## Module shape

The interpolation subsystem becomes one deep Module. Its external Seam consists
of:

- construction with `InterpolationModel` and `BlochOperator`;
- selection of `WignerSeitz()` or `MinimumDistance()` through `real_space`;
- pointwise calculation through `interpolate`;
- path workflows through `band_structure` and, later, `fermi_surface`.

The construction Seam has two real Adapters:

1. the in-memory Wannierization `Model` Adapter;
2. Wannier90 file-data Adapters in `src/io/w90`.

Both must produce the same internal types. File readers must not recreate a
parallel family of tight-binding types.

The expected Depth comes from hiding gauge transformation, quotient Fourier
coefficients, real-space replica selection, symmetry closure, packed evaluation,
derivatives, diagonalization, and observable dependency resolution. This gives
callers Leverage through one learned Interface and gives maintainers Locality by
placing each mathematical rule in one implementation.

## Core data model

### `RealSpaceDomain`

```julia
struct RealSpaceDomain{T}
    vectors::Vector{Vec3{Int}}
    cartesian_vectors::Vector{Vec3{T}}
    vector_index::Dict{Vec3{Int},Int}
end
```

This is the common finite domain on which every primitive operator in one
`InterpolationModel` is represented. It does not remember whether it arose from
Wigner--Seitz or minimum-distance selection. It owns canonical ordering,
fractional/Cartesian consistency, uniqueness, and R-vector-to-index lookup. The
crystal lattice remains owned by `InterpolationModel` rather than being duplicated
here.

All operators are remapped onto the union of the R vectors required by their
real-space scheme, symmetry closure, and Hermitian partners. An individual
operator may therefore contain zero coefficient blocks on part of the common
domain. This shared-domain invariant is preferred because it makes the R axis
uniform and permits one phase matrix for every requested operator. If a future
benchmark demonstrates pathological zero padding, multiple domain groups may be
introduced as a hidden implementation optimization without changing the public
Interface.

### `RealSpaceOperator`

```julia
struct RealSpaceOperator{A,L}
    coefficients::A
    law::L
end
```

The canonical coefficient layout is

```text
n_wannier x n_wannier x component_shape... x n_Rvectors
```

with a dense strided `Array{Complex{T},N}` in the initial implementation. The first
two axes are always Wannier matrix indices, zero or more following axes retain the
logical physical-component shape, and the final axis is the common
`RealSpaceDomain`. Thus a scalar, vector, and rank-two Cartesian tensor use,
respectively,

```text
coefficients[m, n, iR]
coefficients[m, n, alpha, iR]
coefficients[m, n, alpha, beta, iR]
```

The component shape is derived without redundant metadata:

```julia
component_shape(operator) = size(operator.coefficients)[3:(end - 1)]
```

For evaluation, all axes except the last are reshaped without copying to
`(n_wannier^2 * prod(component_shape)) x n_Rvectors` and multiplied by the common
`n_Rvectors x n_kpoints` phase matrix. Outputs are reshaped back to their logical
component axes. A general operator law supplies and validates its component shape;
it is not restricted to Cartesian tensors or powers of three.

Do not store `Vector{Matrix}` or matrices whose entries are `StaticVector`s.
Static vectors remain appropriate for lattice vectors and small symmetry tensors,
not for bulk operator coefficients.

The operator's identity comes from its key in the containing `NamedTuple`; remove
the mutable semantic role currently played by `TBOperator.name::String`.

### `InterpolationModel`

```julia
struct InterpolationModel{C,B,D,O,S}
    crystal::C
    basis::B
    real_space::D
    operators::O
    symmetry::S
end
```

The exact internal types can evolve, but the following invariants are required:

- `operators` is an extensible typed collection, not one field per possible
  physical property;
- `real_space` is the one common `RealSpaceDomain` indexing the final axis of every
  primitive operator;
- `hamiltonian` is one entry in that collection;
- an instance stores only primitive operators supplied or constructed for it;
- derived quantities and temporary buffers are not stored permanently;
- adding an operator changes the collection's value/type, not the definition of
  `InterpolationModel`.

Arrays are heap-backed, so a typed collection of array-owning operator records does
not make the inline struct proportional to the coefficient data. If compilation
time from large `NamedTuple`s becomes measurable, planning may erase types after a
one-time lookup; do not sacrifice the hot evaluation loop to dynamic dispatch
without a benchmark.

## Real-space schemes

Define an unexported abstract supertype and export its concrete values:

```julia
abstract type RealSpaceScheme end
struct WignerSeitz <: RealSpaceScheme end
struct MinimumDistance <: RealSpaceScheme end
```

The internal scheme Interface consumes quotient-lattice coefficients plus crystal
and Wannier-basis metadata and emits weighted infinite-lattice representatives.
Both schemes then enter the same symmetry/Hermiticity closure and packing path.

The minimum-distance implementation must preserve its orbital-pair dependence:
candidate translations and weights are determined from
`R + tau_n - tau_m`. Equal-distance replicas are retained together. The packed
final representation uses the model's common R list, with zeros for matrix
elements that do not use a particular representative; pair-dependent selection
must not leak into the evaluation Interface.

Remove the current staged distinction among `WignerSeitzRspace`, `MDRSRspace`, and
`BareRspace`. There is no public `generate_Rspace` or `simplify` step.

## General-operator construction

### Bloch-space input

Normalize all primitive input to an internal shape

```text
n_bands x n_bands x component_shape... x n_kpoints
```

The physical component axes are preserved for natural indexing and flattened only
inside transformation/evaluation kernels. Diagonal inputs need not be
materialized. The Hamiltonian Adapter views
`model.eigenvalues` as a diagonal scalar operator, declares it Hermitian,
time-reversal even, and scalar under spatial operations, and then invokes exactly
the same gauge-transform/Fourier/scheme/closure pipeline as a user operator.

Construction stages are:

1. validate k-point, band, component, gauge, and law dimensions;
2. transform the operator from the input Bloch basis to the Wannier gauge;
3. compute coefficients on the finite quotient lattice;
4. select weighted infinite-lattice representatives with the real-space scheme;
5. close and project symmetry/Hermitian orbits when symmetry is present;
6. form the union of all required R vectors as one `RealSpaceDomain`;
7. remap and pack every `RealSpaceOperator` onto that common domain.

No stage returns a public intermediate object.

### Transformation laws

Introduce an internal `OperatorLaw` Interface with public, readable constructors:

- `Scalar(time_reversal = Even())`;
- `PolarVector(time_reversal = Even())`;
- `AxialVector(time_reversal = Odd())`;
- a general `CartesianTensor(...)` for advanced operators.

The law supplies the component transformation matrix, the antiunitary conjugation
rule, time-reversal parity, and whether Hermitian conjugation relates opposite
real-space coefficients. Polar/axial determinant factors belong here rather than
in individual observables.

Not every object presently called an operator is a lattice-periodic pointwise
operator. Connections, neighbor-pair matrices, and position-like objects with
inhomogeneous transformation terms require dedicated internal laws or construction
recipes. Do not make the generic law incorrect merely to fit these objects.

## Symmetry architecture

Extract a reusable `WannierSymmetry` from the current localization-focused
`SymmetryConstraint`. It owns only information common to localization and
interpolation:

- space-group operations and antiunitary indicators;
- Wannier centers;
- orbital permutations/mixing matrices;
- orbital-center lattice shifts;
- the data needed to construct Wannier-space transport matrices.

`SymmetryConstraint` should contain or reference `WannierSymmetry` and retain only
the IBZ/FBZ, overlap, Schur-coordinate, and optimization-specific tables required
by localization.

### Symmetry-closed real-space support

Closure acts on coefficient labels

```text
(component_indices..., m, n, R)
```

and not on `R` alone. A transformed lattice vector depends on the spatial action
and on the center shifts of the transported input/output orbitals; a tensor law
may also mix arbitrary tuples of component indices. The symmetry kernel may
linearize that tuple internally, but the stored operator and returned results
retain the logical component axes.

For each primitive operator:

1. seed an orbit with every coefficient produced by the selected real-space
   scheme;
2. apply every unitary and antiunitary symmetry operation;
3. add Hermitian-conjugate partners where required;
4. insert all missing transformed R vectors;
5. project/average coefficients over the complete orbit;
6. merge duplicate labels and pack the final operator.

Minimum-distance ties and their weights are orbit members and must be transported
together. Any later truncation is defined on complete symmetry/Hermitian orbits,
not on individual R vectors or matrix entries.

Keep an internal per-query k-space group projector as an independent reference
Adapter for tests. It is deliberately not the production evaluation path. Tests
compare the reference projector with the one-time real-space projection at random
off-mesh points.

## Evaluation and planning

### Fourier kernel

For the model's common `RealSpaceDomain`, form the phase block once

```math
P_{Rk} = exp(i 2 pi k . R)
```

and evaluate every requested operator's packed coefficients with `mul!`. Process k
points in blocks so the phase matrix and outputs have bounded memory. Operators
with different logical component ranks are reshaped to two-dimensional packed
views without copying. A single-k method uses the same kernel with a one-column
block or an allocation-free matrix-vector specialization.

Derivatives are requested from the same operator rather than stored as separate
`TBHamiltonianGradient` and `TBHamiltonianHessian` objects. Their phase factors are
generated from `(i R_cartesian)` and its products. Compute multiple requested
orders in one traversal or packed multiplication when profiling supports it.

### Observable recipes

Define an internal `Observable` Interface with at least:

```julia
result_name(::BandEnergy) = :band_energy
requirements(::BandEnergy) = (...)
assemble!(output, ::BandEnergy, intermediates, workspace) = ...
```

The concrete observable types are public; dependency tokens and assembly methods
are not. The planner:

1. validates that required primitive operators are present;
2. unions and topologically orders intermediate requirements;
3. constructs one phase block for the common real-space domain;
4. requests only needed derivative orders;
5. schedules exactly one Hamiltonian eigendecomposition per k point;
6. reuses eigenvectors and gauge rotations for every requested observable;
7. allocates the singular-keyed result `NamedTuple`.

Missing data errors must identify both the missing primitive operator and the
observable that requested it.

Initial observable recipes replace the current implementations:

- `BandEnergy`;
- `BandVelocity`;
- `SpinExpectation`, including projection direction as recipe data;
- `BerryCurvature` and its supported formulations;
- `OrbitalMagnetization`, including Fermi energy as recipe data.

Formulation choices belong inside the relevant recipe, not in new
observable-specific model types. For example, two Berry-curvature formulas share
the same public observable concept and declare different internal dependencies.

### Diagonalization

Move batched Hermitian diagonalization behind an internal function and workspace.
Do not add broad `LinearAlgebra.eigen`/`eigen!` methods for vectors of matrices.
Progress reporting is a workflow concern and must not occur inside numerical
kernels.

## Workflows

`band_structure(interpolation_model, kpath)` calls `interpolate` with
`BandEnergy()` and returns a `BandStructure` containing path coordinates, labels,
k points, and band energy. Plotting remains in the Makie extension.

`fermi_surface` and density-of-states calculations can later be migrated onto the
same Interface. They choose a sampling geometry and consume pointwise observable
results; they do not implement Fourier interpolation themselves.

## Source layout

Refactor `src/interpolation` toward the following ownership:

```text
src/interpolation/
    interface.jl                 # exported constructors and interpolate
    types.jl                     # core immutable data types and invariants
    real_space.jl                # schemes, quotient coefficients, packing
    symmetry.jl                  # real-space orbit closure/projection
    fourier_kernel.jl            # packed blocked evaluation
    planning.jl                  # requirements and shared intermediates
    diagonalization.jl           # private batched Hermitian solver
    observables/
        band_energy.jl
        band_velocity.jl
        spin_expectation.jl
        berry_curvature.jl
        orbital_magnetization.jl
    workflows/
        band_structure.jl
```

The exact number of files is secondary to ownership: Fourier evaluation is
implemented once, symmetry closure is implemented once, and observable files
contain formulas rather than storage/evaluation infrastructure.

Update Wannier90 readers to be construction Adapters for `RealSpaceOperator` and
`InterpolationModel`. Readers whose operators contain different R-vector lists
must construct their union and remap them onto one `RealSpaceDomain`. Delete,
rather than wrap, the old public files/types after all callers have moved.

## Implementation sequence

### Phase 0: characterize behavior and performance

- Record numerical reference outputs from current Hamiltonian, spin, position,
  Berry-curvature, and orbital-magnetization tests.
- Add benchmarks for Hamiltonian-only interpolation and a combined
  Hamiltonian/Berry/spin calculation over representative k-point counts.
- Record allocations, phase-construction time, Fourier time, and diagonalization
  time separately.
- Use these references to preserve physics, not old return types or type names.

Verification:

```bash
julia --project=test test/runtests.jl test/interpolation
julia --project=benchmark benchmark/interpolation.jl
```

Create or adjust the benchmark entry point if the second command does not yet
exist.

### Phase 1: core types and ordinary Hamiltonian path

- Rename `KspaceStencil` and `KspaceStencilShells` to `KSpaceStencil` and
  `KSpaceStencilShells`, including constructors, accessors, tests, and docs.
- Implement `RealSpaceDomain`, `RealSpaceOperator`, the two real-space schemes,
  the common-domain invariant, and logical tensor-shaped coefficient storage.
- Implement the general Bloch-to-Wannier gauge transformation and quotient Fourier
  coefficient construction.
- Implement the Hamiltonian input Adapter through the general operator path.
- Implement the packed blocked Fourier kernel.
- Add `InterpolationModel`, `BandEnergy`, `interpolate`, and `band_structure`
  using the new path.
- Test Wigner--Seitz and minimum-distance results against the Phase 0 references.

Exit criterion: ordinary band interpolation and Wannier90 Hamiltonian input work
without constructing any legacy R-space or interpolator type.

### Phase 2: shared planning and derivatives

- Implement observable requirements and the internal execution plan.
- Replace stored Hamiltonian gradient/Hessian operators with derivative requests
  on the Hamiltonian.
- Implement workspace-backed batched diagonalization without extending
  `LinearAlgebra.eigen` for collections.
- Implement `BandVelocity`.
- Verify that requesting `BandEnergy` and `BandVelocity` together performs one
  Fourier evaluation of the Hamiltonian and one eigendecomposition per k point.

Exit criterion: the combined calculation is not slower than two independent
calculations and performs strictly less repeated work.

### Phase 3: symmetry extraction and real-space closure

- Extract `WannierSymmetry` and make localization's `SymmetryConstraint` use it.
- Implement scalar-even Hamiltonian closure first, including antiunitary
  operations and Hermitian partners.
- Extend closure to general component representations and time-reversal parity.
- Integrate minimum-distance equal-replica weights with orbit closure.
- Implement the independent per-query projector used only as a reference.
- Add symmetry diagnostics suitable for construction errors and test failures.

Exit criterion: ordinary Fourier evaluation of closed operators satisfies random
off-mesh covariance and little-group degeneracies to the chosen numerical
tolerance, without per-query group averaging.

### Phase 4: general operators and existing observables

- Implement `BlochOperator` validation and transformation-law constructors.
- Migrate spin and spin projection to `SpinExpectation`.
- Migrate position/connection primitives with their correct dedicated laws.
- Migrate Berry curvature as one observable recipe with formulation data.
- Migrate orbital magnetization and its additional primitive operators.
- Confirm that multi-observable requests share phases, derivatives,
  eigendecompositions, and basis rotations.

Exit criterion: every currently supported physical quantity is accessible through
`interpolate` and no observable-specific interpolator struct remains.

### Phase 5: I/O, workflows, and deletion

- Convert `read_w90_hr`, `read_w90_tb`, `read_w90_spn`, and related readers to the
  new data model.
- Convert writers to consume `RealSpaceOperator` or the relevant operator from an
  `InterpolationModel`.
- Migrate Fermi-energy and Fermi-surface tools to consume `interpolate` results.
- Remove old exports and implementations: `TBOperator`, `AbstractTBInterpolator`,
  `HamiltonianInterpolator`, `SpinInterpolator`,
  `SpinProjectionInterpolator`, `BerryCurvatureInterpolator`,
  `OrbitalMagnetizationInterpolator`, `BareRspace`, `WignerSeitzRspace`,
  `MDRSRspace`, `generate_Rspace`, `simplify`, public `fourier`/`invfourier`, and
  stored Hamiltonian derivative operator types.
- Replace old shallow-module tests with tests at the new Interface. Retain focused
  internal mathematical tests only for Fourier and symmetry kernels.
- Rewrite interpolation documentation and all examples; do not document aliases
  for removed names.

Exit criterion: `rg` finds no use of removed names outside historical material,
all package tests pass, and the documentation builds.

## Correctness tests

Tests at the external Seam must cover:

1. exact reproduction of input-mesh Hamiltonian and general-operator matrices;
2. equality of the Hamiltonian convenience Adapter and an explicitly constructed
   scalar-even `BlochOperator`;
3. Wigner--Seitz and minimum-distance agreement with trusted references;
4. Hermiticity at random on- and off-mesh k points;
5. covariance of every supported transformation-law class at random off-mesh
   points;
6. equality between real-space closure and the independent per-query projector;
7. exact little-group degeneracies for a graphene-like test;
8. antiunitary covariance and Kramers degeneracy for a spinful test;
9. preservation of minimum-distance equal-replica weights under symmetry closure;
10. equality of quantities evaluated alone and in a combined request;
11. clear construction errors for dimension mismatches, missing operator laws, and
    missing primitive dependencies;
12. Wannier90 read/write round trips through the new types.

Use test tolerances derived from input numerical precision and separately report
source-data symmetry noise from interpolation error.

## Performance acceptance criteria

- Steady-state internal `interpolate!` performs no per-k-point heap allocation;
  output arrays and eigensolver workspaces are preallocated.
- All requested operators use one phase matrix for the model's common
  `RealSpaceDomain` per k-point block.
- A combined observable request performs one Hamiltonian eigendecomposition per
  k point.
- Hamiltonian derivatives are generated from R factors and are not duplicated in
  persistent coefficient storage.
- Symmetry-preserving evaluation has no runtime factor proportional to the number
  of symmetry operations; that work is paid once during construction.
- Hamiltonian-only interpolation is no slower than the Phase 0 implementation for
  representative path and dense-grid sizes.
- Combined Berry-curvature/spin/band calculations materially reduce time and
  allocations relative to invoking the old interpolators separately.
- Benchmark both single-threaded and package-standard threaded execution, avoiding
  nested task/BLAS oversubscription.

## Completion checklist

- [x] New public Interface implemented and documented.
- [x] `KSpaceStencil` naming is used consistently in source, tests, and docs.
- [x] Hamiltonian uses the general operator construction path.
- [x] Both real-space schemes emit one common `RealSpaceDomain` for all operators.
- [x] Scalar, vector, matrix, and higher-rank component axes remain logical in
      storage/results and are flattened only in kernels.
- [x] Symmetry-closed R support and coefficient projection implemented.
- [x] Unitary, antiunitary, scalar, polar-vector, axial-vector, and general tensor
      tests pass.
- [x] Observable dependency planning shares all expensive intermediates.
- [x] Dense-grid readiness tests confirm that intermediate storage scales with
      batch size and that observable assembly accepts destination views.
- [x] Existing physical quantities migrated.
- [ ] Wannier90 I/O migrated.
- [ ] Old interpolation Interface deleted with no compatibility wrappers.
- [ ] Focused and full test suites pass.
- [ ] Documentation build passes.
- [ ] Performance acceptance criteria recorded in benchmark output.

Final verification:

```bash
julia --project=test test/runtests.jl
julia --project=docs docs/make.jl
```
