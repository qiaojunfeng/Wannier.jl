# AGENTS.md

## Scope
- These instructions apply to the entire repository.

## Dev environment tips
- Use the latest stable version of Julia.
- From the repo root, instantiate once with:
  - `julia --project -e 'using Pkg; Pkg.instantiate()'`
- Keep local dependencies in sync after dependency changes:
  - `julia --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`
- Run a focused test while iterating on a feature:
  - `julia --project=test test/runtests.jl <test file name>...`
  - Example: `julia --project=test test/runtests.jl w90/win.jl`
- This repo uses workspace projects for `test` and `docs`; run commands with the matching `--project` when needed.

## Testing instructions
- Run package tests from the repo root:
  - `julia --project -e 'using Pkg; Pkg.test()'`
- If you change docs content or APIs, also ensure docs build succeeds:
  - `julia --project=docs docs/make.jl`

## Code and style expectations
- Use 4 spaces for indentation.
- Keep implementations small and type-stable where practical.
- Prefer explicit, readable linear algebra.
- Update docstrings and docs pages when changing public APIs in `src/`.
- Add or update tests for behavior changes, including edge cases for dimensions/units/conventions.

## PR checklist
- Recommended PR title format: `<short summary>`
- Ensure CI-equivalent checks pass locally:
  - Package tests
  - Docs build when docs/public APIs changed
- Keep changes focused; avoid unrelated refactors in the same PR.
- Summarize user-visible API changes in PR description and update README/docs examples when relevant.
- Confirm examples and snippets still run when changing user-facing API behavior.
