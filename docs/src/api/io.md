# Input/Output

This page lists functions for reading and writing files.

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
    * `A` for `amn` matrices
    * `M` for `mmn` matrices
    * `E` for `eig` matrices
    * `S` for `spn` matrices

!!! note

    In most cases, for arrays we adopt the convention that `n_bands` is the first index,
    `n_wann` is the second index, and `n_kpts` is the third index.
    For example, `A` for the `amn` matrices is a 3D array of size `(n_bands, n_wann, n_kpts)`.

```@meta
CurrentModule = Wannier
```

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
read_amn
read_orthonorm_amn
write_amn
read_mmn
write_mmn
read_eig
write_eig
Chk
read_chk
write_chk
get_model
get_A
read_win
read_wout
read_nnkp
write_nnkp
read_unk
write_unk
read_w90
read_w90_post
write_w90
read_w90_band
write_w90_band
get_kpoints
read_w90_wsvec
read_w90_tbdat
read_w90_tb
read_spn
write_spn
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

## Quantum ESPRESSO files

```@docs
read_qe_band
guess_kpath

```
