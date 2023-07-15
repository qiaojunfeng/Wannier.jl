# Foreword to the tutorials

```@meta
CurrentModule = Wannier
```

These tutorials assume you have already been familiar with the
Wannierization workflow using `QE` and `Wannier90`, a good starting
point could be the tutorials of
[`Wannier90`](https://github.com/wannier-developers/wannier90).

The input files of each tutorial can be found in the
[`WannierDocs.jl`](https://github.com/qiaojunfeng/WannierDocs.jl/) repo.
First `git clone` the repo, and run it when you are going through the tutorials.

!!! warning

    In the repo, the convergence parameters (kpoint sampling, ...) in the
    input files (and the generated `amn`/`mmn`/`eig`/... files) are extremely
    loose! This is on purpose to reduce storage size and save computational time.
    Thus, it is expected that the band interpolation, the resolution of real space
    WFs are not good.

    You should rerun the calculations with production-quality convergence
    parameters, and the results should be much better.
