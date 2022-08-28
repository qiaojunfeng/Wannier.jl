# Wannier.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qiaojunfeng.github.io/Wannier.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://qiaojunfeng.github.io/Wannier.jl/dev)
[![CI](https://github.com/qiaojunfeng/Wannier.jl/workflows/CI/badge.svg)](https://github.com/qiaojunfeng/Wannier.jl/actions?query=workflow%3ACI)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![codecov](https://codecov.io/gh/qiaojunfeng/Wannier.jl/branch/main/graph/badge.svg?token=J2c9HRdk59)](https://codecov.io/gh/qiaojunfeng/Wannier.jl)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

**A playground for experimentation with Wannier functions (WFs)[^MV97][^MMYSV12].**

## Features

### Wannierization

- maximal localization for isolated bands, e.g. insulators
  - different from[^MV97], optimize on unitary matrix manifolds (adaptation of [^DLL19] to isolated bands)
- disentanglement for entangled bands, e.g. metal
  - different from[^SMV01], optimize on Stiefel manifolds[^DLL19]
- parallel transport gauge[^GLS19]
  - you can further improve the spread by optimization w.r.t. a single rotation matrix[^QMP21]
- split valence and conduction WFs from a valence + conduction calculation[^QMP21]
  - as a by-product, automated initial projection for valence or conduction WFs
  - for the initial projection of valence + conduction calculation, you can start with either
    conventional spdf projection, SCDM[^DL18], or an automated projection and disentanglement
    from pseudopotential orbitals[^QPM21]
  - different from SCDM, the valence+conduction manifold is chosen by the valence+conduction calculation,
    instead of SCDM μ and σ. Moreover, works in reciprocal space thus more memory-efficient
- constraint center for max localization or disentanglement[^QMP21]
  - similar to[^WLPMM14], add an Lagrange multiplier term to spread functional, but optimize
    on matrix manifolds, and applying to both max localization and disentanglement
    (whereas in [^WLPMM14] the center is constrained only during max localization)

### Interpolation

Two algorithms:

- Wigner-Seitz (WS) interpolation
- Minimal-distance replica selection (MDRS) method

for band structure along a kpath or on a grid.

### Real space WFs

- output `xsf` or `cube` file
- evaluate operators in real space

## Installation

Install with the Julia package manager [Pkg](https://pkgdocs.julialang.org/),
just like any other registered Julia package:

```jl
pkg> add Wannier  # Press ']' to enter the Pkg REPL mode.
```

or

```jl
julia> using Pkg; Pkg.add("Wannier")
```

### CLI

Additionally, there is a command line interface

```bash
julia --project deps/build.jl install  # install CLI
```

The executable will be installed in ```~/.julia/bin/wannier```.
After appending `~/.julia/bin` to your `$PATH`, you can use the CLI as follows:

```bash
$ wannier -h


  wannier v0.1.0

Julia package for Wannier functions.

Usage

  wannier <command>
...
```

Note since julia needs precompilation of the code, there will be some delay when running the CLI.

## Usage

For example, running a maximal localization can be easily achieved by

```jl
using Wannier

model = read_w90("silicon")
A = max_localize(model)
write_amn("silicon.amn", A)
```

For a much more detailed overview, please see
[the User Guide documentation](https://qiaojunfeng.github.io/Wannier.jl/stable/user/).

## Development

- clone repo: `git clone https://github.com/qiaojunfeng/Wannier.jl`
- install pre-commit: `pre-commit install`
- test:

  ```bash
  julia --project=.  # start REPL
  ]                  # activate Pkg mode
  test               # run tests
  ```

## Contributing

The code initially started with Antoine Levitt's repo
[wannier](https://github.com/antoine-levitt/wannier), and went through a series of
refactorization, bug fixes, and feature additions.

This is a research code mainly for development and testing.
Issues and pull requests are welcome!

## References

[^MV97]: Marzari, N. & Vanderbilt, D.,
    Maximally localized generalized Wannier functions for composite energy bands,
    [Physical Review B, 1997, 56, 12847](https://doi.org/10.1103/physrevb.56.12847)
[^MMYSV12]: Marzari, N.; Mostofi, A. A.; Yates, J. R.; Souza, I. & Vanderbilt, D.,
    Maximally localized Wannier functions: Theory and applications,
    [Reviews of Modern Physics, 2012, 84, 1419](https://doi.org/10.1103/revmodphys.84.1419)
[^SMV01]: Souza, I.; Marzari, N. & Vanderbilt, D.,
    Maximally localized Wannier functions for entangled energy bands,
    [Physical Review B, 2001, 65, 035109](https://doi.org/10.1103/physrevb.65.035109)
[^DLL19]: Damle, A.; Levitt, A. & Lin, L.,
    Variational Formulation for Wannier Functions with Entangled Band Structure,
    [Multiscale Modeling & Simulation, 2019, 17, 167](https://doi.org/10.1137/18m1167164)
[^GLS19]: Gontier, D.; Levitt, A. & Siraj-dine, S.,
    Numerical construction of Wannier functions through homotopy,
    [Journal of Mathematical Physics, 2019, 60, 031901](https://doi.org/10.1063/1.5085753)
[^QPM21]: Qiao, J.; Pizzi, G. & Marzari, N.,
    Projectability disentanglement for accurate high-throughput Wannierization,
    xxx
[^QMP21]: Qiao, J.; Marzari, N. & Pizzi, G.,
    Automated separate Wannierization for valence and conduction manifolds,
    xxx
[^DL18]: Damle, A. & Lin, L.,
    Disentanglement via Entanglement: A Unified Method for Wannier Localization
    [Multiscale Modeling & Simulation, 2018, 16, 1392](https://doi.org/10.1137/17m1129696)
[^WLPMM14]: Wang, R.; Lazar, E. A.; Park, H.; Millis, A. J. & Marianetti, C. A.,
    Selectively localized Wannier functions,
    [Physical Review B, 2014, 90, 165125](https://doi.org/10.1103/physrevb.90.165125)
