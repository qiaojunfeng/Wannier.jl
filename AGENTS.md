# AGENTS.md

## Scope
- These instructions apply to the entire repository.
- `CLAUDE.md` is a symlink to this file; edit only `AGENTS.md`.

## Development
- Use the latest stable Julia release.
- From the repo root, instantiate once with `julia --project -e 'using Pkg; Pkg.instantiate()'`.
- After dependency changes, keep the environment in sync with `julia --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`.
- While iterating, run focused tests with `julia --project=test test/runtests.jl <test file name>...`.
  - Example: `julia --project=test test/runtests.jl test/symmetry/operations.jl`.
- This repo uses separate projects for `test` and `docs`; use the matching `--project` flag when needed.

## Testing
- Run the full package test suite from the repo root with `julia --project=test test/runtests.jl`.
- If you change docs content or public APIs, also verify the docs build with `julia --project=docs docs/make.jl`.
- In tests, call non-exported APIs with module qualification, for example `Wannier.project_covariant!(...)`.

## Naming
Name verbosity scales with distance from context: the closer a name sits to
the user-facing API, the more of its meaning it must carry itself; the deeper
it sits inside a mathematical kernel, the shorter it should be, so that the
structure of the mathematics stays visible. Three tiers:

1. **Exported API** (functions, types, keyword arguments): full words, no
   abbreviations, no Unicode — users must be able to guess, type, and grep
   them. Examples: `symmetry_constraint`, `SymmetrizedModel`,
   `clean_littlegroup_reps!`, `atol_degeneracy`. Functions are verbs (with
   `!` for mutation), types are nouns.
2. **Internal structure** (unexported helpers, struct fields, cross-file
   plumbing): concise domain abbreviations are welcome when they are either
   community-standard (`mmn`, `amn`, `kpb`, `ibz`, `fbz`, `nnkp`) or defined
   in the glossary below / the file-header comment of the module that owns
   them (`sc`, `ws`, `Lmat`, `ikb`). An abbreviation that is neither is a
   bug: spell it out or add it to the glossary.
3. **Mathematical kernels** (function bodies, hot loops, derivation-heavy
   code): very short names, including Unicode, are preferred when they
   mirror a citable derivation — the paper/SM or the documents under
   `unfold/` — so the code reads like the equations (`Ω`, `ε`, `q`, `𝒦`,
   `M̃` as `Mt`). Anchor rule: the surrounding comment or docstring must
   point to the equation or document section that fixes the notation.
   Unicode never appears in exported names or file names.

Bridge rule: when one object crosses tiers, record the correspondence once,
in the header comment of the owning file (e.g. "`ws.Mt` = M̃, the
Wannier-gauge overlap of SM Eq. (S21) = what the user-facing docs call the
projected overlap matrix"). Renames of public names require prior approval.

### Glossary of blessed abbreviations
- `sc`  — `SymmetryConstraint`
- `ws`  — workspace struct of the current kernel
- `sb`  — `SchurBasis`
- `ibz` / `fbz` — irreducible / full Brillouin zone
- `kpb` — k+b neighbor bookkeeping (w90 convention)
- `Mt`, `Dt` — M̃ (Wannier-gauge overlaps), D̃ (phase-dressed orbital rep)
- `Lmat`, `Rmat` — the 𝓛/𝓡 orbital transports of the SM
- w90 file stems (`mmn`, `amn`, `eig`, `nnkp`, `chk`, `spn`, `uHu`) and
  their IBZ variants (`immn`, `iamn`, `ieig`, `isym`)

## Code Style
- Use 4 spaces for indentation.
- Keep implementations small and type-stable where practical.
- Prefer explicit, readable linear algebra.
- Update docstrings and docs pages when changing public APIs in `src/`.
- Add or update tests for behavior changes, including edge cases for dimensions/units/conventions.
- If a local variable has the same name as a keyword argument, Julia lets you omit the keyword name in the call, for example `foo(x; y)`.

## PR checklist
- Recommended PR title format: `<short summary>`
