# Input/Output

```@meta
CurrentModule = Wannier
```

The reading and writing functions are implemented in the
[`WannierIO.jl`](https://github.com/qiaojunfeng/WannierIO.jl)
package. However, here are also some convenience functions which
wrap the corresponding functions in `WannierIO.jl`, to utilize the
`struct`s defined in `Wannier90.jl`, e.g. [`KspaceStencil`](@ref),
[`RGrid`](@ref), etc.

!!! tip

    In most cases, the units of the function arguments and returns are in
    angstrom unit for lattice, and fractional w.r.t lattice for
    atomic positions, etc.

!!! tip

    The following abbreviations are used throughout the code and documentation:
    * `n_bands` for number of bands
    * `n_wann` for number of WFs
    * `n_kpts` for number of kpoints
    * `n_bvecs` for number of b-vectors
    * `n_atoms` for number of atoms
    * `U` for `amn` or the gauge matrices
    * `M` for `mmn` matrices
    * `E` for `eig` matrices
    * `S` for `spn` matrices

!!! note

    In most cases, for arrays we adopt the convention that `n_bands` is the first index,
    `n_wann` is the second index, and `n_kpts` is the third index.
    For example, `U` for the gauge matrices is a 3D array of size `(n_bands, n_wann, n_kpts)`.

## Contents

```@contents
Pages = ["io.md"]
Depth = 2
```

## Index

```@index
Pages = ["io.md"]
```

## Wannier90 files

```@autodocs
Modules = [Wannier]
Pages   = [
    "io/w90/amn.jl",
    "io/w90/band.jl",
    "io/w90/chk.jl",
    "io/w90/model.jl",
    "io/w90/nnkp.jl",
    "io/w90/tb.jl",
]
```

## File manipulation

### Truncate Wannier90 matrices

!!! tip

    Here are some functions to remove some bands from `mmn`, `eig`, or `UNK` files,
    so as to skip rerunning NSCF calculations and `pw2wannier90.x`.

```@autodocs
Modules = [Wannier]
Pages   = ["io/truncate.jl"]
```

## 3D visualization files

```@autodocs
Modules = [Wannier]
Pages   = [
    "io/volume/xsf.jl",
    "io/volume/bxsf.jl",
    "io/volume/cube.jl",
]
```

## Interface to DFT codes

```@autodocs
Modules = [Wannier]
Pages   = [
    "io/interface/mud.jl",
]
```
