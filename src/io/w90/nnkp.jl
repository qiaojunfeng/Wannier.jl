export read_nnkp_compute_bweights, write_nnkp

"""
    $(SIGNATURES)

Read the `nnkp` file.

This function calls `WannierIO.read_nnkp` to parse the file, compute the bweights
of b-vectors, and returns a [`KspaceStencil`](@ref) (while `WannierIO.read_nnkp` only
returns a `NamedTuple`).
"""
function read_nnkp_compute_bweights(filename::AbstractString)
    nnkp = WannierIO.read_nnkp(filename)
    return KspaceStencil(
        nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G
    )
end

"""
    $(SIGNATURES)

Write nnkp that can be used by `pw2wannier90`.

# Arguments
- `filename`: the filename to write to
- `kstencil`: a [`KspaceStencil`](@ref) object

!!! tip

    Some important tags in `nnkp` file (can be passed as keyword arguments):
    - `n_wann`: the number of WFs, needed by `pw2wannier90`
    - `exclude_bands`: the bands (often semicore states) to exclude, needed by
        `pw2wannier90`

    For other keyword arguments, see [`WannierIO.write_nnkp`](@ref).
"""
function write_nnkp(filename::AbstractString, kstencil::KspaceStencil; kwargs...)
    return WannierIO.write_nnkp(
        filename;
        lattice=real_lattice(reciprocal_lattice(kstencil)),
        recip_lattice=reciprocal_lattice(kstencil),
        kstencil.kpoints,
        kstencil.kpb_k,
        kstencil.kpb_G,
        kwargs...,
    )
end
