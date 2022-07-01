# Wannier.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://qiaojunfeng.github.io/Wannier.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://qiaojunfeng.github.io/Wannier.jl/dev)
[![CI](https://github.com/qiaojunfeng/Wannier.jl/workflows/CI/badge.svg)](https://github.com/qiaojunfeng/Wannier.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/qiaojunfeng/Wannier.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/qiaojunfeng/Wannier.jl)

**A playground for experimentation with Wannier functions.**

## Features

* Wannierization
  * maximal localization for isolated bands, e.g. insulators
  * disentanglement for entangled bands, e.g. metal
  * parallel transport gauge
  * split valence and conduction WFs from a valence + conduction calculation
* interpolation
  * band structure

## Installation

Install with the Julia package manager [Pkg](https://pkgdocs.julialang.org/), just like any other registered Julia package:

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

model = read_seedname("silicon")
A = max_localize(model)
write_amn("silicon.amn", A)
```

---

For a much more detailed overview, please see [the User Guide documentation](https://qiaojunfeng.github.io/Wannier.jl/stable/user/).

## Contributing

This is a research code mainly for development and testing.
Issues and pull requests are welcome!
