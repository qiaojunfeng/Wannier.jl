# Benchmark performance

## Usage

From the repository root, instantiate the benchmark environment once:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
```

Run the complete package benchmark suite with the reproducible single-threaded
configuration used by `PkgBenchmark`:

```bash
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    julia --project=benchmark benchmark/runbenchmarks.jl
```

The interpolation group can be run directly, either single-threaded or in a
representative four-thread configuration while keeping BLAS single-threaded to
avoid nested oversubscription:

```bash
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=benchmark -e \
    'using BenchmarkTools; include("benchmark/benchmarks.jl"); run(SUITE["interpolation"])'

JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 julia --project=benchmark -e \
    'using BenchmarkTools; include("benchmark/benchmarks.jl"); run(SUITE["interpolation"])'
```

`SUITE["interpolation"]["band energy"]` contains source-mesh, 256-point path,
and ``32^3`` dense-grid Hamiltonian-only workloads for both correctness and
performance tracking. The `band-energy stages` subgroup separates phase
construction, the Fourier sum, and diagonalization. The remaining groups
compare combined observable requests with independent calls.

Inspired by <https://github.com/JuliaFolds/Transducers.jl/tree/master/benchmark>.
