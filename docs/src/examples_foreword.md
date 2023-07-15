# Foreword to the examples

```@meta
CurrentModule = Wannier
```

If you have already been familiar with the Wannierization workflow using
Wannierization workflow using `Quantum ESPRESSO` and `Wannier90`, the following
examples should be easy to follow.
As a side note, some good [`Wannier90`] tutorials can be found at

- examples folder of [`Wannier90`](https://github.com/wannier-developers/wannier90/tree/develop/examples)
- <https://github.com/wannier-developers/wannier-tutorials>

!!! warning

    In this repo, the convergence parameters (kpoint sampling, ...) when generating
    the input files (e.g., `amn`/`mmn`/`eig`/...) are quite loose!
    This is on purpose to reduce storage size and save computational time.
    Thus, it is expected that the band interpolation might have some defects.

    In real applications, one should always run calculations with production-quality
    convergence parameters for high-quality band interpolations.

## WannierDatasets

Throughout the examples, we need to load input files to construct `Model`s for
Wannierizations.
Often, preparing these input files requires density-functional-theory
calculations which can be time-consuming, and distracting ourselves from
the main purpose of Wannierizations and Wannier interpolations.

Julia has a convenient [`artifact`](https://pkgdocs.julialang.org/v1/artifacts/)
system allowing us to load binaries or data files easily in Julia script or REPL,
without the need of manually downloading and placing them in the right folder.

In `Wannier.jl`, we use such system to load our pre-computed
[`WannierDatasets`](https://github.com/qiaojunfeng/WannierDatasets), which
contains some typical materials, e.g., silicon, copper, graphene, etc.
This allows easy loading of input files and construct Wannier functions hassle-free.

The `Wannier.Datasets` submodule provides several functions to help inspect and
load the datasets.
First, we need to `using` the submodule to bring the functions into scope

```@repl intro_artifacts
using Wannier.Datasets
```

List all available datasets by

```@repl intro_artifacts
list_datasets()
```

Load a dataset into a `Wannier.jl` [`Model`](@ref) by

```@repl intro_artifacts
model = load_dataset("Si2")
```

The dataset is just a folder containing input files, you can inspect the folder by

```@repl intro_artifacts
show_dataset("Si2")
```

Finally, you can directly access the individual files by `artifact` string

```@repl intro_artifacts
artifact"Si2/si2.win"
readlines(artifact"Si2/si2.win")
```
