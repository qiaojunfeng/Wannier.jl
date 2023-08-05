"""
    $(SIGNATURES)

Compute
``\\langle u_{m, k + b_1}| H | u_{n, k + b_2} \\rangle`` matrices
from overlap matrices
``\\langle u_{m, k} | u_{n, k + b} \\rangle``
and Hamiltonian.

Use the following formula:

```
\\langle u_{m, k + b_1}| H | u_{n, k + b_2} \\rangle
= \\sum_{i,j} \\langle u_{m, k + b_1} | u_{i, k} \\rangle
    \\langle u_{i, k} | H | u_{j, k} \\rangle
    \\langle u_{j, k} | u_{n, k + b_2} \\rangle
= M_{k, b}^{\\dagger} H_{k} M_{k, b}
```

"""
function compute_uHu(
    overlaps::AbstractVector,
    kpb_k::AbstractVector,
    kpb_G::AbstractVector,
    eigenvalues::AbstractVector,
)
    WannierIO._check_dimensions_M_kpb(overlaps, kpb_k, kpb_G)
    nkpts = length(kpb_k)
    @assert nkpts == length(eigenvalues) "wrong n_kpoints in eigenvalues"

    uHu = map(zip(overlaps, eigenvalues)) do (Mₖ, εₖ)
        # note I need to use reshape to convert it into a row vector,
        # cannot use transpose since that one is recursive
        (adjoint.(Mₖ) .* Ref(Diagonal(εₖ))) * reshape(Mₖ, (1, :))
    end
    return uHu
end
