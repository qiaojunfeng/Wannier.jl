# Getting Started

## Installation

### Julia

First install [Julia](https://julialang.org/), you can,

- download official release from <https://julialang.org/downloads/>
- or use [`juliaup`](https://github.com/JuliaLang/juliaup) which can better handle multiple Julia versions

!!! tip

    Quick notes for installing `juliaup`, on Linux or macOS:

    ```bash
    curl -fsSL https://install.julialang.org | sh
    juliaup add release
    ```

### Wannier.jl

#### From Julia package manager

Install with the Julia package manager [Pkg](https://pkgdocs.julialang.org/),
just like any other registered Julia package:

```julia
pkg> add Wannier  # Press ']' to enter the Pkg REPL mode.
```

or

```julia
julia> using Pkg; Pkg.add("Wannier")
```

This is the recommended way to install `Wannier.jl` for end users.

#### From GitHub

```bash
git clone https://github.com/qiaojunfeng/Wannier.jl.git
cd Wannier.jl
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

#### From tarball

```bash
tar xvf Wannier.jl.tgz
cd Wannier.jl
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

!!! tip

    If you install with `Pkg`, then start Julia REPL just with `julia`;
    if you install with git or tarball, then start with `julia --project=PATH_OF_Wannier.jl`.

### Command-line interface

Additionally, there is a command-line interface (CLI)

```bash
cd Wannier.jl
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

!!! note

    Since julia needs precompilation of the code, there will be some delay when running the CLI.

## Usage

For example, running a maximal localization can be easily achieved by

```julia
using Wannier

model = read_w90("silicon")
A = max_localize(model)
write_amn("silicon.amn", A)
```

For a much more detailed usage, please refer to the API documentation
and the source code of each function.

## Development

The repo is hosted at <https://github.com/qiaojunfeng/Wannier.jl>

### Pre-commit

To ensure uniform code style, we use [pre-commit](https://pre-commit.com/),
which auto check and fix code style before each `git commit`. Install the hooks with

```bash
cd Wannier.jl
pre-commit install
```

### Tests

Run tests with

```bash
cd Wannier.jl
julia --project=.  # start REPL
]                  # activate Pkg mode
test               # run tests
```
