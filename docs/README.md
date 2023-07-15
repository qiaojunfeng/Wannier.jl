# Documentation

## Writing

The documentation is written in Markdown, then processed by `Documenter.jl`.

### Examples

1. The examples are written as Julia scripts, in folder `examples/`
2. The scripts are processed by `Literate.jl`, converted to
    - Markdown
    - Jupyter notebook
    - Julia script without comments
3. The Markdowns are processed by `Documenter.jl` to generate HTMLs

The quickest starting point is looking at existing examples,
e.g. [`maximal_localization.jl`](src/examples/maximal_localization.jl),
or [`band_structure.jl`](src/examples/band_structure.jl) for integration with `PlotlyJS`.

## Build

See [`make_serve.sh`](./make_serve.sh).
