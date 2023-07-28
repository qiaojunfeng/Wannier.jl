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
    kpoints = nnkp.kpoints
    recip_lattice = nnkp.recip_lattice
    kpb_k = nnkp.kpb_k
    kpb_G = nnkp.kpb_G
    n_bvecs = length(kpb_k[1])

    # Generate bvectors from 1st kpoint, in fractional coordinates
    bvectors = zeros(Vec3{Float64}, n_bvecs)
    ik = 1
    for ib in 1:n_bvecs
        ikpb = kpb_k[ik][ib]
        G = kpb_G[ik][ib]
        bvectors[ib] = recip_lattice * (kpoints[ikpb] + G - kpoints[ik])
    end

    bweights = compute_bweights(bvectors)
    kgrid_size = guess_kgrid_size(kpoints)
    return KspaceStencil(
        recip_lattice, kgrid_size, kpoints, bvectors, bweights, kpb_k, kpb_G
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
