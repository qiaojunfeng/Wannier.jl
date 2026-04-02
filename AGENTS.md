# AGENTS.md

## Scope
- These instructions apply to the entire repository.

## Development
- Use the latest stable Julia release.
- From the repo root, instantiate once with `julia --project -e 'using Pkg; Pkg.instantiate()'`.
- After dependency changes, keep the environment in sync with `julia --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`.
- While iterating, run focused tests with `julia --project=test test/runtests.jl <test file name>...`.
  - Example: `julia --project=test test/runtests.jl w90/win.jl`.
- This repo uses separate projects for `test` and `docs`; use the matching `--project` flag when needed.

## Testing
- Run the full package test suite from the repo root with `julia --project -e 'using Pkg; Pkg.test()'`.
- If you change docs content or public APIs, also verify the docs build with `julia --project=docs docs/make.jl`.
- In tests, call non-exported APIs with module qualification, for example `WannierIO.write_HH_R(...)`.

## Code Style
- Use 4 spaces for indentation.
- Keep implementations small and type-stable where practical.
- Prefer explicit, readable linear algebra.
- Update docstrings and docs pages when changing public APIs in `src/`.
- Add or update tests for behavior changes, including edge cases for dimensions/units/conventions.
- If a local variable has the same name as a keyword argument, Julia lets you omit the keyword name in the call, for example `foo(x; y)`.

## PR checklist
- Recommended PR title format: `<short summary>`
- Ensure CI-equivalent checks pass locally:
  - Package tests
  - Docs build when docs/public APIs changed
- Keep changes focused; avoid unrelated refactors in the same PR.
- Summarize user-visible API changes in PR description and update README/docs examples when relevant.
- Confirm examples and snippets still run when changing user-facing API behavior.
