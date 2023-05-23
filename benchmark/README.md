# Benchmark performance

## Usage

Start julia REPL in current directory `julia --project`, then `instantiate` the
project if it is the first time you run the benchmark,

```julia
using Pkg; Pkg.instantiate()
```

Afterwards, run

```julia
include("runbenchmarks.jl")
```

Inspired by <https://github.com/JuliaFolds/Transducers.jl/tree/master/benchmark>.
