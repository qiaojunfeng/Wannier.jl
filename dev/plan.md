## Plan: Modernize Localization Architecture (REVISED)

Modernize Wannier.jl by separating optimization concerns from model data, introducing a `Problem`-first architecture, and preparing packed array storage for CPU and future GPU while prioritizing low-risk refactors first. **Key refinement (user decision):** Remove `Cache` struct from one-off calls; keep Cache only in optimization closures where it amortizes allocation (~1,400 buffers over 200 iterations). Use array-type dispatch (`Array`, future `CuArray`) as the default device distinction. Fix MU mutation corruption early (Phase 2).

**Steps**
1. Phase 1: Establish architectural baseline and scope boundaries. ← **NEXT (Restarted)**
   - Confirm that `Model` remains physics/input-only and does not store optimization terms or layout.
   - Re-establish baseline with focused localization tests and full package test command.
   - Re-apply shared MU/UtMU kernel extraction in a fresh Commit B, preserving exact behavior.
   - Re-check load-order discipline for new modules.

2. Phase 2: Layout Utilities & Fix Gradient Corruption.
   - **Commit C**: Extract layout conversions from disentangle.jl to `src/common/layouts.jl`:
     - X_Y_to_U!, U_to_X_Y, XY_to_X_Y!, X_Y_to_XY, GU_to_GX_GY
     - Frozen-band masking utilities
   - **Commit D**: Remove Cache from one-off calls + fix MU corruption:
     - Remove `Cache` instantiation from `omega(bvectors, M, U)`, `omega_grad(penalty, bvectors, M, U)`, `center(bvectors, M, U)`
     - Inline buffer allocation (MU, UtMU, G, r) in these functions; call kernels directly
     - Fix `omega_grad!` side-effect that mutates `cache.MU` — use temporary scratch matrix, keep MU read-only
     - Cache now exists ONLY in optimization closures (disentangle, max_localize, coopt) where it matters
     - Validation: All 19 tests pass, convergence identical

3. Phase 3: Array-Dispatch Consolidation (Preparatory).
   - **Commit E**: Standardize array-type dispatch contracts (no behavior change):
     - Keep kernels and helper interfaces generic over array container type
     - Ensure CPU path uses packed arrays as primary internal representation
     - Validation: Tests pass, no numerical drift

4. Phase 4: Introduce LocalizationProblem Struct (*depends on phases 2–3*).
   - **Commit F**: Create `LocalizationProblem` and unified objective builder:
     - Struct: model, penalty_terms, solver_options
     - **NOTE: Problem does NOT own Cache**; optimization closures/solver path allocate workspace directly from array types
     - Implement `build_fg!(problem, U_or_XY)` objective builder
     - Keep existing entry points working; internally construct Problem + call solver
     - Validation: Same 19 tests pass, convergence identical

5. Phase 5: Migrate Optimization Entry Points (*depends on phase 4*).
   - **Commit G**: Migrate max_localize to Problem-based:
     - Create `localize_isolated_bands(model, penalty, solver="Optim")` (new API)
     - Old `max_localize` becomes thin wrapper (backward-compat fallback)
   - **Commit H**: Migrate disentangle equivalently
   - **Commit I**: Migrate magnetic/coopt variants
   - Validation: New API produces identical results to old

6. Phase 6: CPU Packed-Storage Finalization (*GPU deferred*).
   - **Commit J**: Finalize packed CPU storage (4D/3D arrays):
     - Extend `compute_MU_UtMU!` dispatch for packed arrays (Vector dispatch already exists)
     - Add accessor/reshape helpers
     - Benchmark: allocation count, cache locality vs nested containers
     - If packed layout wins, make it default and remove unnecessary nested-layout paths
     - Validation: Benchmarks confirm benefit

**Execution order and commit strategy**
1. **Commit A** (Phase 1): Recreate baseline test checkpoint from clean restart state
2. **Commit B** (Phase 1): Re-extract shared MU/UtMU kernel in fresh commit (behavior-preserving, load-order verified)
3. **Commit C** (Phase 2): Extract layout conversions & frozen utilities from disentangle.jl → src/common/layouts.jl
4. **Commit D** (Phase 2): Remove Cache from one-off calls, inline buffer allocation; fix MU corruption in omega_grad! (keep cache only in opt closures)
5. **Commit E** (Phase 3): Standardize array-type dispatch contracts (no layout/backend trait layer)
6. **Commit F** (Phase 4): Introduce LocalizationProblem struct + build_fg! objective builder; workspace determined by array types in solver/closure path (NOT Problem)
7. **Commit G** (Phase 5): Migrate max_localize → localize_isolated_bands (new API, old is wrapper)
8. **Commit H** (Phase 5): Migrate disentangle equivalently
9. **Commit I** (Phase 5): Migrate magnetic/coopt variants
10. **Commit J** (Phase 6): Packed CPU storage finalization (array-type dispatch kernels, benchmarks)

**Key principle**: Each commit is focused, test-gated, behavior-preserving. Plan is restarted from Phase 1; treat all commits A through J as pending until revalidated in sequence.

**Relevant files**
- `/home/jqiao/git/Wannier.jl/src/model.jl` — keep `Model` data-only; do not embed optimization layout/terms.
- `/home/jqiao/git/Wannier.jl/src/spread.jl` — extract shared kernels and normalize objective/gradient interfaces.
- `/home/jqiao/git/Wannier.jl/src/localization/max_localize.jl` — migrate to Problem-based objective construction.
- `/home/jqiao/git/Wannier.jl/src/localization/disentangle.jl` — migrate layout conversions and `fg!` generation into shared abstractions.
- `/home/jqiao/git/Wannier.jl/src/localization/coopt.jl` — migrate magnetic co-optimization objective/gradient wiring.
- `/home/jqiao/git/Wannier.jl/src/localization/constrain_center/coopt.jl` — migrate constrained-center magnetic objective wiring.
- `/home/jqiao/git/Wannier.jl/src/localization/opt_rotate.jl` — align any inline objective/gradient computation with shared pipeline.
- `/home/jqiao/git/Wannier.jl/src/Wannier.jl` — include new modules for kernels/layout helpers/problem abstractions.
- `/home/jqiao/git/WannierIO.jl/src/w90/mmn.jl` — refactor `Mmn` storage to packed arrays or add packed representation API.
- `/home/jqiao/git/WannierIO.jl/src/w90/hr_dat.jl` — evaluate packed operator storage (`H`) with clear accessor helpers.
- `/home/jqiao/git/WannierIO.jl/src/w90/tb_dat.jl` — evaluate packed storage for `H` and position operators (`r_x`, `r_y`, `r_z`).
- `/home/jqiao/git/Wannier.jl/src/io/w90/model.jl` — adapt model construction path from parser output to packed-internal layout with minimal conversion overhead.

**Verification**
1. Run full test suite from repo root: `julia --project -e 'using Pkg; Pkg.test()'`.
2. Add/execute focused localization tests: `julia --project=test test/runtests.jl localization/disentangle.jl localization/max_localize.jl` (or existing equivalent files).
3. Add finite-difference checks for each term-composition path (standard spread, center constraint, spin coupling).
4. Add round-trip tests for each layout: U -> XY -> U and gradient mapping consistency.
5. Add benchmark comparisons for old nested storage vs packed storage in CPU kernels.

**Decisions**
- ✅ Use a dedicated Problem struct; no convenience wrapper required.
- ✅ Keep optimization layout (U/XY/spin-XY) inside Problem, not Model.
- ✅ Backward compatibility is not a constraint; prioritize coherent architecture.
- ✅ **Cache struct role**: Keep only in optimization closures (amortizes ~1,400 allocations over 200 iterations). Remove from one-off omega/omega_grad/center calls (inline allocation). Buffer strategy is determined by concrete array types, not a backend/layout tag in Problem.
- ✅ **Fix MU corruption early** (Phase 2, Commit D): omega_grad! currently mutates cache.MU as side-effect; refactor to use scratch matrix and keep MU read-only.
- ✅ CPU path should move toward packed multidimensional storage where dimensions are uniform.
- ✅ Share CPU/GPU math kernels through one array-based kernel layer and array-type-specific allocation/adaptation only.
- ✅ Preserve mathematical readability by adding named accessor/helper APIs over packed storage.
- ✅ Refactor WannierIO after Wannier.jl stabilizes (separate, later PR).
- ✅ GPU readiness: rely on array-type dispatch (`Array` now, `CuArray` later) with GPU kernels deferred.
- ✅ Public API naming redesigned for clarity; remove legacy names directly (no compatibility wrappers).
- ✅ Include solver adapter placeholder (thin abstraction layer); Optim.jl default, alternatives (Manopt/Optimization) deferred after numerical parity validation.

**Further Considerations**
1. **Load-order discipline**: kernels.jl, layouts.jl, and future problem.jl must have clear include order. Document dependencies in src/Wannier.jl comment.
2. **Closure pattern template**: disentangle.jl closure (reusing Cache across iterations) is reference pattern for max_localize/coopt migration in Phase 5.
3. **Solver buffer strategy**: Different solvers may prefer different buffer policies (stack vs async queue vs GPU async), but the default differentiation remains array types (`Array` / `CuArray`) unless policy control becomes a concrete need.