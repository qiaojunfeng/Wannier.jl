"""
    $(SIGNATURES)

Compute the so-called `uHu` matrix, LVTS12 Eq. 96:
``\\langle u_{m, k + b_1}| H | u_{n, k + b_2} \\rangle``,
approximately from the overlap matrices
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

This is an approximation to the true `uHu` matrix that should be computed
e.g. in the DFT plane-wave basis.

Note the `overlaps` and `eigenvalues` should be the outputs from a DFT to
Wannier interface code, i.e., they should contain `n_bands`-sized matrices
instead of Wannier-gauge `n_wannier`-sized matrices, to minimize the effect
of band truncation.
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

"""
Compute H_kb = U'_k H_k M_kb U_kb, LVTS12 Eq. 86 and 91.

```
\\langle u_{m, k}| H | u_{n, k + b} \\rangle
= \\sum_{i} \\langle u_{m, k} | H | u_{i, k} \\rangle
    \\langle u_{i, k} | u_{n, k + b} \\rangle
= H_{k} M_{k, b}
```

In principle, this should be computed in the DFT plane-wave basis, similar to
[`compute_uHu`](@ref); however, no DFT codes implement this, so we approximate
it using the Hamiltonian and overlap matrices.

Note the `overlaps` and `eigenvalues` should be the outputs from a DFT to
Wannier interface code, i.e., they should contain `n_bands`-sized matrices
instead of Wannier-gauge `n_wannier`-sized matrices, to minimize the effect
of band truncation.
"""
function compute_hamiltonian_times_position_kspace(
    overlaps::AbstractVector,
    kpb_k::AbstractVector,
    kpb_G::AbstractVector,
    eigenvalues::AbstractVector,
)
    WannierIO._check_dimensions_M_kpb(overlaps, kpb_k, kpb_G)
    nkpts = length(kpb_k)
    @assert nkpts == length(eigenvalues) "wrong n_kpoints in eigenvalues"

    Hkb = map(zip(eigenvalues, overlaps)) do (εₖ, Mₖ)
        Ref(Diagonal(εₖ)) .* Mₖ
    end
    # Hkb can be indexed by Hkb[ik][ib][m, n]
    return Hkb
end
