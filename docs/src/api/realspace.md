# Real space WFs

This page lists functions for processing WFs in real space.

Normally operators are computed in reciprocal space, but sometimes it might be useful to
evaluate them in real space. For example, computing higher moment of WFs.

```@meta
CurrentModule = Wannier
```

## Contents

```@contents
Pages = ["realspace.md"]
Depth = 2
```

## Index

```@index
Pages = ["realspace.md"]
```

## Real-space grid struct

```@autodocs
Modules = [Wannier]
Pages   = ["common/rgrid.jl"]
```

## Read/write real-space WFs

```@autodocs
Modules = [Wannier]
Pages   = ["realspace/wavefunction.jl"]
```

## Evaluate operators in real space

```@autodocs
Modules = [Wannier]
Pages   = ["realspace/moment.jl"]
```
