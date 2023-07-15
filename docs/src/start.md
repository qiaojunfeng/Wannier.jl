# Getting Started

## Installation

### Julia

First install [Julia](https://julialang.org/), you can,

- use [`juliaup`](https://github.com/JuliaLang/juliaup)
- or download releases from <https://julialang.org/downloads/>

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
pkg> add Wannier  # Press ']' to enter the Pkg mode
```

or

```julia
julia> using Pkg; Pkg.add("Wannier")
```

#### From GitHub

If you want to play with the code, you can clone the repo:

```bash
git clone https://github.com/qiaojunfeng/Wannier.jl.git
cd Wannier.jl
julia --project=. -e 'using Pkg; Pkg.update(); Pkg.resolve(); Pkg.instantiate()'
```

#### From tarball

Download releases from <https://github.com/qiaojunfeng/Wannier.jl/releases>

```bash
tar xvf Wannier.jl.tar.gz
cd Wannier.jl
julia --project=. -e 'using Pkg; Pkg.update(); Pkg.resolve(); Pkg.instantiate()'
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

The executable will be installed in ```~/.julia/bin/wjl```.

After appending `~/.julia/bin` to your `$PATH`, you can use the CLI as follows:

```bash
$ wjl -h


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
U = max_localize(model)
write_amn("silicon.amn", U)
```

For a much more detailed usage, please refer to the API documentation
and the source code of each function.

## Development

The repo is hosted at <https://github.com/qiaojunfeng/Wannier.jl>

To reduce the startup latency of the `Wannier.jl` package,
and to allow smoother user/development experience, some of the
functionality are implemented inside standalone packages:

- The input/output functions are inside
    [`WannierIO.jl`](https://github.com/qiaojunfeng/WannierIO.jl) repo,
    and its documentation at
    [io.wannierjl.org](https://io.wannierjl.org/)
- The plotting related code is inside
    [`WannierPlots.jl`](https://github.com/qiaojunfeng/WannierPlots.jl) repo,
    and its documentation at
    [plots.wannierjl.org](https://plots.wannierjl.org/)

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
