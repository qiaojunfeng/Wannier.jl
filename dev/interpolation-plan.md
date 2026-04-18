# Interpolation API Simplification + Performance Plan

Companion to [dev/plan.md](plan.md) (localization rewrite). Scope: `src/interpolation/*` only.

## Priorities

1. **API cleanliness** — one canonical arg order, one entry point per transform direction, no overload sprawl.
2. **Performance** — precomputed phase matrix, BLAS3 core, thread over independent axis, share phase across interpolators.
3. **Collapse boilerplate** — unify `TBHamiltonianGradient` / `TBHamiltonianHessian` / etc. as R-weighted views; trivial wrapper interpolators deleted.

Backward compatibility not a constraint beyond a one-release deprecation shim.

## Problems found

### API inconsistency

- `fourier(kpts, op, Rspace)` vs `invfourier(Rspace, op, kpts)` — arg order flipped.
  [src/interpolation/fourier.jl:47-85](../src/interpolation/fourier.jl#L47), [src/interpolation/fourier.jl:130-178](../src/interpolation/fourier.jl#L130).
- Callback variants also flipped: `f(ik, iR, phase)` vs `f(iR, ik, phase)`.
  [src/interpolation/fourier.jl:202-239](../src/interpolation/fourier.jl#L202).
- Four entry points per direction: `AbstractVector{<:AbstractMatrix}`, `Array{T,3}`, `TBOperator`, callback.

### Hot-loop perf smells

[src/interpolation/fourier.jl:59-66](../src/interpolation/fourier.jl#L59) and [:142-150](../src/interpolation/fourier.jl#L142):

- `exp(±im·2π·k·R)` recomputed every (k,R). N_k · N_R exp calls per operator, repeated per interpolator.
- `Oᴿ .+= phase * Oᵏ` — elementwise broadcast, no BLAS `axpy!` / `gemm`.
- `fill!` inside outer loop (minor).
- No threading; embarrassingly parallel over R (k→R) and k (R→k).

### Upstream redundancy (larger impact than fourier itself)

- Velocity FD: 6× `invfourier!` same R-op, shifted k.
  [src/interpolation/hamiltonian_gradient.jl:161-162](../src/interpolation/hamiltonian_gradient.jl#L161).
- Effective-mass FD: 19× `invfourier` same op.
  [src/interpolation/hamiltonian_hessian.jl:189-190](../src/interpolation/hamiltonian_hessian.jl#L189).
- BerryCurvature: 4× invfourier, same kpts, different ops.
  [src/interpolation/berry_curvature.jl:104](../src/interpolation/berry_curvature.jl#L104) and :128, :155, :189.
- `compute_D_matrix` re-invfouriers H already computed by caller.
  [src/interpolation/position.jl:95](../src/interpolation/position.jl#L95), [src/interpolation/position.jl:234](../src/interpolation/position.jl#L234).
- Trivial wrapper interpolators: [src/interpolation/hamiltonian_gradient.jl:20-21](../src/interpolation/hamiltonian_gradient.jl#L20).

## Key design decisions

1. **Unified arg order both directions**: `(Rspace, operator, kpoints)`.
2. **Single entry point per direction**: `AbstractVector{<:AbstractMatrix}` is canonical; drop `Array{T,3}` overloads (provide reshape helper) and callback variants (replace by `phase_matrix`).
3. **Precomputed phase matrix**: `Φ[iR, ik] = cis(±2π · R·k)` as `Matrix{ComplexF64}`. One `cis` per (R,k), reusable across every interpolator hitting the same kpts.
4. **BLAS3 core**: operator vector reshaped to `N_wann² × N_R` complex matrix; transform is one `gemm`.
5. **KPointCache**: opt-in struct carrying `Φ`, eigenvalues, gauge `U` across interpolators, and across FD stencil offsets.
6. **Lazy R-weighted wrapper**: `TBRWeighted(op, α)` replaces `TBHamiltonianGradient` / `TBHamiltonianHessian` etc. — no new storage.

## Plan

### Phase 1 — Unify API (breaking)

1. Standardize both directions: `fourier(Rspace, operator_k, kpoints)` / `invfourier(Rspace, operator_R, kpoints)`.
2. Drop `Array{T,3}` overloads; document `reshape_operator` helper.
3. Remove callback variants; replace with `phase_matrix(Rspace, kpts) -> Matrix{ComplexF64}`.
4. Clean up `TBOperator` shorthand — single `invfourier(tb::TBOperator, kpts)` definition.
5. Deprecation shims for one release.

### Phase 2 — BLAS3 core

1. Precompute `Φ[iR, ik] = cis(±2π · R·k)` once per call.
2. Reshape operator vector → `N_wann² × N_R` dense complex matrix. Transform:
   - k→R: `Or_flat = Ok_flat * Φ' / N_k`
   - R→k: `Ok_flat = Or_flat * Φ`
   Single `gemm`, single allocation.
3. Scalar fallback path for non-Matrix eltypes (Vec3 for Berry Ω).
4. Land as 5-line no-API-change patch first (zero risk) — prove speedup before breaking API.
5. Expected: 10–50× speedup for N_wann ≥ 20.

### Phase 3 — Thread + share phase

1. `@threads` on scalar fallback outer axis.
2. `KPointCache` struct: `Φ`, eigenvalues, gauges `U` for a set of kpts. Passed into every interpolator so Position / Berry / Spin / Mag reuse a single factorization.
3. FD stencils (velocity, effective mass) build `Φ` once over all stencil kpts, not per offset.

### Phase 4 — Collapse boilerplate

1. `TBHamiltonianGradient`, `TBHamiltonianHessian`, ... are all `R^α ⊗ op_R`, α ∈ {0,1,2}. Unify as lazy `TBRWeighted(op, α)`.
2. `compute_D_matrix` accepts pre-computed H(k) + eigen; stop re-invfouriering
   [src/interpolation/position.jl:234-236](../src/interpolation/position.jl#L234).
3. Delete trivial `Interpolator(H) = Interpolator(H, TB...(H))` wrappers; make interpolators callable on raw `TBOperator`.

### Phase 5 — Tests + bench

1. Keep existing tests green via deprecation shims (one release).
2. Bench harness vs `postw90` for H(k), v(k), Ω(k) on silicon; commit baseline numbers.
3. Add phase-matrix equivalence tests (new BLAS path vs legacy scalar loop, tolerance `atol=1e-12`).

## Sequencing / risk

- **Phase 2 first** as zero-risk patch on current API — prove speedup.
- **Phases 1 + 2 together** as single breaking release.
- **Phases 3, 4** additive thereafter.
- **Phase 5** continuous.

## Risk log

| Risk | Phase | Mitigation |
|------|-------|-----------|
| Breaking external users of `fourier` / `invfourier` | 1 | One-release deprecation shim; coordinate with WannierIO.jl consumers |
| BLAS path only strided dense | 2 | Fallback scalar path for Vec3-valued / sparse operators |
| `KPointCache` adds hidden state | 3 | Opt-in; plain calls still work |
| Threading races with shared `Φ` | 3 | `Φ` is read-only after construction |
| FD stencil refactor alters velocity numerics | 3, 4 | Parity gate vs current velocity/effective-mass reference values |

## Out-of-scope

- `src/interpolation/Rspace.jl` type definitions stay as-is.
- `src/interpolation/operator.jl` `TBOperator` struct unchanged.
- `src/interpolation/fermi_energy.jl`, `berry_curvature.jl` physics kernels — only their invfourier call sites touched.
- Localization rewrite — see [dev/plan.md](plan.md).
