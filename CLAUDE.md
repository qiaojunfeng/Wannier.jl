# CLAUDE.md

## Scope
- These instructions apply to the entire repository.

## Development
- Use the latest stable Julia release.
- From the repo root, instantiate once with `julia --project -e 'using Pkg; Pkg.instantiate()'`.
- After dependency changes, keep the environment in sync with `julia --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`.
- While iterating, run focused tests with `julia --project=test test/runtests.jl <test file name>...`.
  - Example: `julia --project=test test/runtests.jl test/test1.jl`.
- This repo uses separate projects for `test` and `docs`; use the matching `--project` flag when needed.

## Testing
- Run the full package test suite from the repo root with `julia --project=test test/runtests.jl`.
- If you change docs content or public APIs, also verify the docs build with `julia --project=docs docs/make.jl`.
- In tests, call non-exported APIs with module qualification, for example `MyPackage.read_special_file(...)`.
