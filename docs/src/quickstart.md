# Quick start

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

For a much more detailed usage, please refer to the API documentation
and the source code of each function.

## Development

- clone repo: `git clone https://github.com/qiaojunfeng/Wannier.jl`
- install pre-commit: `pre-commit install`
- test:

  ```bash
  julia --project=.  # start REPL
  ]                  # activate Pkg mode
  test               # run tests
  ```
