# Wannier.jl

A Julia package for playing with Wannier functions (WFs), including

- Generation of WFs from density functional theory (DFT) calculations
- Wannier interpolation of operators, e.g. band structure

The purpose of this package is to write clean and modularized functions for Wannierization
and Wannier interpolation, together with the expressive and interactive Julia language,
to enable fast prototyping and development of new ideas.

Compared with [Wannier90](http://www.wannier.org/), the package provides different
Wannierization algorithms using matrix manifold optimization.
Currently, MPI is not supported, so probably not suitable for large systems.

The package is still under development, and the API is subject to change.

## Features

- Wannierization

  - maximal localization for isolated bands, e.g. insulators
  - disentanglement for entangled bands, e.g. metal
  - parallel transport gauge
  - split valence and conduction WFs
    - automated initial projection for valence or conduction WFs
  - constrain WF center

- Interpolation of operators, e.g. band structure
- Real space WFs

  - output `xsf` or `cube` file
  - evaluate operators in real space
