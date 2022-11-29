# Input/Output

```@meta
CurrentModule = Wannier
```

The reading and writing functions are implemented in the
[`WannierIO.jl`](https://github.com/qiaojunfeng/WannierIO.jl)
package. However, here are also some convenience functions which
wrap the corresponding functions in `WannierIO.jl`, to utilize the
`struct`s defined in `Wannier90.jl`, e.g. [`BVectors`](@ref),
[`RGrid`](@ref), etc.

!!! tip

    In most cases, the units of the function arguments and returns are in angstrom unit for lattice,
    and fractional w.r.t lattice for atomic positions, etc.

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

```@docs
read_orthonorm_amn
read_w90_band
write_w90_band
write_chk
read_w90
read_w90_interp
write_w90
read_nnkp
write_nnkp
read_w90_tb
```

## File manipulation

### Truncate Wannier90 matrices

!!! tip

    Here are some functions to remove some bands from `mmn`, `eig`, or `UNK` files,
    so as to skip rerunning NSCF calculations and `pw2wannier90.x`.

```@docs
truncate_mmn_eig
truncate_unk
truncate_w90
```

## 3D visualization files

```@docs
read_xsf
write_xsf
read_cube
write_cube
```

## Misc

```@docs
KPathInterpolant
get_symm_idx_label
Model(chk::WannierIO.Chk)
```
