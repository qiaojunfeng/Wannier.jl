# Why do we need this folder?

For some functions, e.g., `read_nnkp`, they are already defined in `WannierIO.jl`;
however, in `Wannier.jl`, we have a wrapper function (`read_nnkp_compute_bweights`),
which return a `KspaceStencil` struct instead of a `NamedTuple`, also containing the
weights of bvectors. This makes the code more user-friendly.

Thus, the functions in the `io/w90` folder only contains wrapper functions
which accept or return structs defined in this package.
All the actual parsing codes are strictly in `WannierIO.jl`.
