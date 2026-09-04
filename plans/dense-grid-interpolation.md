# Dense-grid Wannier interpolation

## Status and dependency

The core redesign is complete; its current interface is documented in the
[interpolation API](../docs/src/api/interpolation.md) and
[Fourier-interpolation theory](../docs/src/theory/fourier.md). This follow-on
adds regular-grid, streaming, and FFT capabilities to the implemented
`InterpolationModel`, common `RealSpaceDomain`, observable recipes, and
batch-oriented primitive/assembly flow.

The core redesign deliberately keeps dense-grid concerns out of its public
Interface. It must, however, satisfy the dense-grid readiness contract in that
plan. This plan then adds the public and internal behavior needed for regular
meshes containing approximately one million or more k points.

## Goals

1. Interpolate observables on large regular k meshes with bounded temporary
   memory.
2. Avoid materializing k-point coordinates, Fourier phases, Hamiltonians,
   derivatives, eigensystems, or observable intermediates for the complete mesh.
3. Support collected in-memory output, caller-owned or memory-mapped output, and
   streaming Brillouin-zone reductions through the same calculation path.
4. Exploit symmetry by evaluating irreducible representatives when the requested
   observable has a well-defined, degeneracy-safe transformation law.
5. Provide an optimized blocked direct transform for arbitrary point sets and an
   FFT implementation for compatible regular grids.
6. Select batch size and Fourier implementation from a memory/cost model rather
   than requiring routine manual tuning.
7. Share phases, primitive operator interpolation, derivatives, eigensystems, and
   basis rotations across all requested observables in each batch.
8. Make working memory scale with batch size rather than total k-point count.

## Non-goals

- Change `InterpolationModel`, `RealSpaceDomain`, `RealSpaceOperator`, or the
  observable recipe Interface established by the core redesign.
- Implement a general distributed-memory or multi-node execution system.
- Implement out-of-core multidimensional FFTs in the first version.
- Guarantee symmetry reduction for quantities that are not individually defined
  at degeneracies.
- Add GPU execution in the first version. Packed storage and batch kernels must
  remain compatible with a future accelerator Adapter.
- Require callers to choose a Fourier implementation or batch size in ordinary
  use.

## Public Interface

### Lazy regular grid

```julia
grid = KPointGrid(
    (100, 100, 100);
    shift = (0, 0, 0),
)
```

`KPointGrid` is an indexable, lazy collection of fractional k-point coordinates.
For dimensions `(N1, N2, N3)`, zero-based grid indices `j`, and a shift expressed
in grid-step units, it represents

```math
k(j) = ((j_1 + s_1)/N_1, (j_2 + s_2)/N_2, (j_3 + s_3)/N_3).
```

The type stores dimensions, shift, length, and a documented linear ordering. It
does not allocate a `3 x n_kpoints` coordinate array. It satisfies the same
indexing/range-access Interface already accepted by `interpolate`, so arbitrary
vectors of k points and lazy grids use the same external Seam.

Required accessors are full-word exported names:

```julia
grid_dimensions(grid)
grid_shift(grid)
```

The dense query grid is distinct from `KSpaceStencil`: `KSpaceStencil` describes
the source Wannierization mesh and its finite-difference neighbors, whereas
`KPointGrid` is a requested evaluation geometry.

### Collected result

```julia
result = interpolate(
    model,
    grid,
    BerryCurvature();
    memory_limit = 1_000_000_000,
)

result.berry_curvature
```

`memory_limit` bounds temporary working memory and excludes the returned result.
Before allocating, `interpolate` estimates the result size. If it is unusually
large relative to available memory, the error/warning must recommend
`interpolate!` with caller-owned or memory-mapped storage rather than failing
after partial work.

The default memory limit is conservative and documented in bytes. The planner
derives batch size from this limit; routine users do not choose a batch size.

### Caller-owned or memory-mapped result

Promote `interpolate!` from the internal core-redesign function to the public
Interface:

```julia
destination = (
    berry_curvature = Array{Float64}(undef, 3, n_wannier, length(grid)),
)

interpolate!(
    destination,
    model,
    grid,
    BerryCurvature();
    memory_limit = 1_000_000_000,
)
```

The destination is a `NamedTuple` whose singular keys, element types, component
axes, band axes, and final k-point axis must match the observable result schema.
Any suitable `AbstractArray` is accepted, including a memory-mapped array. The
method validates all destinations before beginning evaluation and returns the
destination.

Do not expose output-sink implementation types merely to support arrays. The
public `interpolate!` method is the array Adapter at the output Seam.

### Streaming Brillouin-zone reductions

Calculations that need only an integral must not construct the pointwise field:

```julia
integral = brillouin_zone_integral(
    model,
    grid,
    BerryCurvature(occupation = occupation),
)
```

`brillouin_zone_integral` consumes observable batches immediately, applies grid or
irreducible-point weights, and stores only per-task accumulators. Observable-
specific workflows such as anomalous Hall conductivity may call this function and
add physical prefactors and units; they must not implement interpolation loops.

The reduction convention, Brillouin-zone normalization, units, occupation model,
and treatment of partially occupied degeneracies must be documented explicitly.

## Batch execution Module

### Data flow

Every implementation follows the same flow:

```text
k-point source range
    -> primitive real-space operator evaluation
    -> shared derivatives and eigensystem
    -> observable assembly
    -> output sink or reduction accumulator
```

Only the first and last steps vary between collected output, caller-owned output,
and reductions. Observable formulas consume a common internal batch
representation and do not know whether primitive operators came from a direct
transform or an FFT.

### Memory model

For a direct-transform batch of size `B`, the principal complex working arrays
scale approximately as

```math
16 N_R B
+ 16 N_W^2 N_{primitive} B
+ M_{eigensolver}
+ M_{observable},
```

where `N_primitive` includes physical components and requested derivative orders.
The planner computes an exact byte estimate from array element types and logical
component shapes. It reserves space for:

- the phase block;
- packed primitive outputs;
- derivative-weighted outputs;
- eigensolver matrices and workspaces;
- observable intermediates;
- one or two batch buffers when pipelining is enabled.

The selected batch size must:

- fit within `memory_limit` including alignment and task-local workspaces;
- be at least one, otherwise report the minimum required memory;
- favor matrix dimensions efficient for BLAS without exceeding the limit;
- remain independent of the total grid size;
- be capped where larger blocks no longer improve measured throughput.

An internal batch-size override is permitted for tests and benchmarks but is not
part of the initial public Interface.

### Output sinks

Define one internal sink Interface only because there are multiple real Adapters:

- collected array sink used by `interpolate`;
- caller-owned array sink used by `interpolate!`;
- reduction sink used by `brillouin_zone_integral`.

A sink accepts a result schema, a global k-point range, and one completed
observable batch. Sink code must not trigger Fourier evaluation or
diagonalization. A callback/HDF5 sink can be added later only when a concrete
second external-storage use case requires it.

### Workspaces and reentrancy

Plans contain immutable dependencies, layouts, and cost decisions. Mutable
scratch arrays live in explicit task-local workspaces. Do not put reusable mutable
scratch storage inside `InterpolationModel`.

The same model can therefore be used concurrently by independent interpolation
calls. Within one call, each active task owns its eigensolver and observable
workspace. Output ranges are disjoint.

## Direct Fourier implementation

### Arbitrary point sets

The blocked direct implementation remains the universal fallback:

1. generate one phase matrix for the common `RealSpaceDomain` and current k-point
   batch;
2. reshape each required operator/derivative coefficient block to two dimensions
   without copying;
3. multiply packed coefficients by the phase matrix with `mul!`;
4. restore logical tensor component axes through views;
5. assemble and write the observable batch.

Where profitable, concatenate required coefficient rows into one planned packed
block so that one large multiplication replaces several small multiplications.
Packing happens once per plan, not once per k-point batch.

### Regular-grid phase recurrence

For `KPointGrid`, avoid evaluating a complex exponential independently for every
`(R,k)` pair. Generate phases in the documented grid order using multiplicative
increments. To control accumulated roundoff:

- compute the exact starting phase for every batch;
- use recurrence only within the batch;
- periodically reseed long inner recurrences when required by accuracy tests.

The recurrence implementation must agree with direct `exp` phase construction at
the selected numerical tolerance for shifted and unshifted grids, negative R
vectors, and R vectors outside the source-grid representative cell.

## FFT implementation

### Eligibility

The FFT implementation is eligible only when:

- the query is a complete `KPointGrid`;
- its ordering and shift are supported;
- all primitive results needed simultaneously fit the memory limit;
- the planner's cost estimate predicts a benefit over direct evaluation;
- requested symmetry reduction does not make direct evaluation on the much
  smaller irreducible set cheaper.

Otherwise the planner silently selects the blocked direct implementation. A
qualified diagnostic function may report the decision and memory estimate for
benchmarking; backend-choice types are not initially exported.

### Mapping arbitrary R vectors

For grid dimensions `N`, each integer R vector maps to a periodic FFT index
`mod(R, N)`. Several minimum-distance or symmetry-generated R vectors can collide
at one FFT index and their coefficient contributions must be summed. A shifted
query grid contributes the corresponding shift phase before the transform.

This mapping makes no assumption that `RealSpaceDomain` is a Wigner--Seitz cell,
rectangular grid, or subset of the source-grid representatives.

For derivatives, multiply coefficients by the actual Cartesian R components
before mapping modulo the dense query grid. Colliding modes generally have
different R factors, so derivative-weighted coefficients must be accumulated
separately; they cannot be recovered by differentiating a prematurely aliased
coefficient grid.

### Storage and execution

Flatten Wannier and logical component axes into transform channels. Perform
batched multidimensional FFTs over channel blocks, then restore component axes for
observable assembly.

The first implementation uses FFT only when the complete set of simultaneously
required primitive k-space arrays fits in working memory. Do not silently spill
large FFT intermediates to disk. Out-of-core FFTs or slab/pencil decompositions
require a separate plan if profiling later justifies their complexity.

Direct and FFT implementations must produce the same primitive-batch convention,
normalization, signs, grid ordering, shift convention, and derivative units.

## Symmetry-reduced dense grids

### Irreducible mapping

When `InterpolationModel` contains symmetry, construct an internal irreducible
mapping for a compatible `KPointGrid`:

- validate that every symmetry operation maps the shifted grid onto itself;
- select deterministic representatives;
- record full-to-irreducible indices and the operation used for reconstruction;
- record multiplicity/integration weights;
- include antiunitary operations where mathematically applicable.

All mapping checks should use integer grid arithmetic after validating the shift,
rather than repeated floating-point nearest-neighbor searches.

### Observable eligibility

Each observable recipe declares whether its result can be reconstructed from an
irreducible representative and provides the required transformation. The planner
uses symmetry reduction only when every jointly requested observable is eligible.

Examples:

- band energy is invariant and eligible;
- total or occupied-subspace Berry curvature has an axial-vector,
  time-reversal-odd transformation and is eligible;
- an individual band's Berry curvature at a degeneracy is not uniquely defined
  and is not automatically eligible;
- general tensor results use their declared `OperatorLaw`/observable law.

If eligibility is absent, evaluate the full grid rather than guessing a law. Keep
an internal switch that forces full-grid evaluation for verification.

### Reconstruction and integration

For a full destination, evaluate representatives in batches and transform each
result directly into its full-grid destination indices. Do not first allocate a
second complete irreducible result unless required by ordering.

For reductions, evaluate each representative once and multiply by its exact grid
multiplicity/weight. Do not reconstruct the full field. Use deterministic pairwise
or compensated summation for numerically sensitive integrals.

The planner compares full-grid FFT cost with irreducible direct-transform cost;
symmetry reduction is not assumed to be faster in every regime.

## Parallel execution

Use one controlled parallel layer at a time:

- direct Fourier multiplication may use multithreaded BLAS for sufficiently large
  blocks;
- small Hermitian eigensystems and observable assembly may use Julia tasks with
  one workspace per task;
- FFT execution uses an explicitly configured FFT thread count;
- never combine task-level batch parallelism with uncontrolled multithreaded BLAS
  or FFT oversubscription.

Benchmark at least two strategies on representative systems:

1. sequential batches with threaded BLAS/FFT, followed by threaded eigensystems;
2. concurrent batch tasks with single-threaded BLAS and independent workspaces.

Select from measured dimensions and costs rather than the total CPU thread count
alone. Results and reductions must be deterministic within documented floating-
point tolerances.

## Source ownership

Extend the interpolation implementation without moving observable formulas:

```text
src/interpolation/
    kpoint_grid.jl              # lazy KPointGrid and grid arithmetic
    batching.jl                 # memory estimates and batch ranges
    sinks.jl                    # array and reduction Adapters
    symmetry_sampling.jl        # irreducible grid mapping/reconstruction
    fourier_kernel.jl           # blocked direct implementation
    fourier_fft.jl              # regular-grid FFT implementation
    integration.jl              # generic Brillouin-zone reduction
```

If the core redesign chooses different filenames, preserve the same ownership:
grid geometry, memory planning, output sinks, Fourier implementations, and
symmetry sampling each have one home. No observable file should contain its own
dense-grid traversal.

## Implementation sequence

### Phase 0: baseline and cost model inputs

- Select representative models spanning small and moderate Wannier dimensions,
  R-domain sizes, scalar/vector/tensor operators, and unitary/antiunitary
  symmetry.
- Benchmark the core blocked implementation at `10^3`, `10^4`, `10^5`, and
  `10^6` arbitrary k points.
- Measure phase generation, packed multiplication, diagonalization, observable
  assembly, output writing, allocations, and peak resident memory separately.
- Record final-output sizes for band energy, total Berry curvature, and
  band-resolved Berry curvature.

Exit criterion: the memory/cost model is based on measured array sizes and timing
data rather than fixed guessed batch sizes.

### Phase 1: lazy grid and bounded direct execution

- Implement `KPointGrid` with allocation-free scalar indexing and efficient range
  materialization into a caller-provided batch buffer.
- Enforce the core plan's batch-oriented primitive/assembly flow.
- Implement exact workspace byte estimation and automatic batch selection.
- Add the public `memory_limit` keyword.
- Verify that direct interpolation over one million points has bounded temporary
  memory and no per-point allocation.

Exit criterion: a million-point Berry-curvature calculation can run without
materializing full-mesh intermediates; only the selected destination scales with
the complete grid.

### Phase 2: destinations and reductions

- Promote and document public `interpolate!`.
- Validate named destination schemas before computation.
- Test ordinary arrays and memory-mapped arrays.
- Implement the internal sink Interface and `brillouin_zone_integral` reduction.
- Add deterministic weighted accumulation and observable-specific workflow hooks.

Exit criterion: pointwise output can be streamed to caller-owned storage and a
Brillouin-zone integral uses memory independent of total grid size.

### Phase 3: symmetry reduction

- Implement exact shifted-grid symmetry compatibility checks and irreducible
  mappings.
- Add observable eligibility and reconstruction methods.
- Implement direct reconstruction into full destinations and weighted reduction
  without reconstruction.
- Verify unitary, antiunitary, scalar, axial-vector, and general-tensor cases.
- Benchmark symmetry reduction against full-grid direct evaluation.

Exit criterion: eligible observables agree with full-grid evaluation while
performing primitive interpolation only at irreducible representatives.

### Phase 4: optimized regular-grid direct transform

- Implement batch-local phase recurrence with controlled reseeding.
- Pack jointly required primitive/derivative coefficient rows once per plan.
- Tune batch dimensions and parallel strategy from representative benchmarks.
- Retain exact-exponential phase generation as an independent test reference.

Exit criterion: recurrence and packing improve throughput without changing results
beyond documented floating-point tolerances.

### Phase 5: FFT implementation and selection

- Add the FFT dependency or extension and reusable FFT plans.
- Implement modulo-grid coefficient accumulation, collisions, shifts, and
  derivative-weighted grids.
- Implement channel batching and logical tensor-axis restoration.
- Add a memory/cost-based implementation selector.
- Compare FFT, direct regular-grid, and symmetry-reduced direct results and timing.

Exit criterion: FFT is selected only for eligible cases where it fits the memory
limit and provides a measured benefit; fallback requires no caller action.

### Phase 6: workflows, documentation, and tuning

- Migrate dense Fermi-surface, density-of-states, Berry-curvature integration, and
  anomalous-Hall workflows to the shared streaming path.
- Document output sizes and when to choose `interpolate!` or an integral workflow.
- Document thread control and reproducible benchmark commands.
- Add performance regression baselines without encoding machine-specific absolute
  timing as correctness tests.

Exit criterion: examples cover a collected path, a million-point memory-mapped
Berry-curvature field, and a symmetry-reduced streaming integral.

## Correctness tests

Tests through the external Seam must cover:

1. `KPointGrid` coordinates, length, ordering, shift, and range access;
2. equality between lazy-grid and explicitly collected k-point inputs;
3. equality across many forced internal batch sizes, including one and the full
   small test grid;
4. equality of `interpolate` and `interpolate!` for scalar, vector, rank-two, and
   higher-rank component shapes;
5. memory-mapped destination equality and complete overwrite of every output
   element;
6. direct recurrence agreement with exact exponential phases;
7. FFT/direct agreement for shifted and unshifted grids;
8. FFT collisions from R vectors differing by dense-grid periods;
9. correct derivatives when colliding modes have different actual R vectors;
10. symmetry-reconstructed full output versus direct full-grid output;
11. irreducible weighted reductions versus full-grid reductions;
12. antiunitary and tensor reconstruction;
13. safe fallback for observables that are not degeneracy-safe;
14. deterministic repeated reductions within documented tolerances;
15. validation failures before partial writes for insufficient memory, malformed
    destinations, incompatible shifted grids, and unsupported result laws.

Internal tests may force direct, recurrence, FFT, full-grid, and irreducible
implementations. These controls remain qualified and unexported.

## Memory and performance acceptance criteria

- Temporary memory is `O(batch_size)` for the direct implementation and does not
  grow with total k-point count.
- `interpolate!` adds no complete-result allocation beyond caller-owned storage.
- Streaming reductions use `O(number_of_tasks)` accumulators and no pointwise
  complete-grid result.
- Steady-state batch kernels perform no per-k-point heap allocation.
- One phase block is shared by every requested primitive operator in a direct
  batch.
- One Hamiltonian eigendecomposition is performed per evaluated k point, not per
  observable.
- Symmetry-reduced integrations interpolate only irreducible representatives.
- FFT plans and work buffers are reused across channel batches and calls within
  one execution plan.
- The planner never selects an implementation whose estimated temporary memory
  exceeds `memory_limit`.
- Million-point benchmarks report total runtime, evaluated points per second,
  peak resident memory, allocated bytes, and time fractions for phases/FFT,
  primitive interpolation, diagonalization, observable assembly, and output.
- Performance reports distinguish total-grid points from actually evaluated
  irreducible points.

## Completion checklist

- [x] Core interpolation redesign and dense-grid readiness contract complete.
- [ ] Lazy `KPointGrid` implemented and documented.
- [ ] Automatic memory estimation and batch sizing implemented.
- [ ] Public `interpolate!` supports ordinary and memory-mapped arrays.
- [ ] Streaming Brillouin-zone reduction implemented.
- [ ] Symmetry-reduced evaluation and reconstruction implemented for eligible
      observables.
- [ ] Regular-grid recurrence implementation verified.
- [ ] FFT implementation handles arbitrary R vectors, collisions, shifts, and
      derivatives.
- [ ] Automatic implementation selection verified against forced references.
- [ ] Scalar, vector, matrix, and higher-rank observable shapes tested.
- [ ] Million-point memory and performance benchmarks recorded.
- [ ] Dense workflows and documentation migrated.
- [ ] Focused and full test suites pass.
- [ ] Documentation build passes.

Final verification:

```bash
julia --project=test test/runtests.jl test/interpolation
julia --project=test test/runtests.jl
julia --project=docs docs/make.jl
```
