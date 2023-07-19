export read_nnkp_compute_weights, write_nnkp

"""
    $(SIGNATURES)

Read the `nnkp` file.

This function calls `WannierIO.read_nnkp` to parse the file, compute the weights
of b-vectors, and returns a [`KgridStencil`](@ref) (while `WannierIO.read_nnkp` only
returns a `NamedTuple`).
"""
function read_nnkp_compute_weights(filename::AbstractString)
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

    weights = compute_weights(bvectors)
    kgrid_size = guess_kgrid_size(kpoints)
    kgrid = KpointGrid(recip_lattice, kgrid_size, kpoints)
    return KgridStencil(kgrid, bvectors, weights, kpb_k, kpb_G)
end

"""
    $(SIGNATURES)

Write nnkp that can be used by `pw2wannier90`.

# Arguments
- `filename`: the filename to write to
- `kstencil`: a [`KgridStencil`](@ref) object

!!! tip

    Some important tags in `nnkp` file (can be passed as keyword arguments):
    - `n_wann`: the number of WFs, needed by `pw2wannier90`
    - `exclude_bands`: the bands (often semicore states) to exclude, needed by
        `pw2wannier90`

    For other keyword arguments, see [`WannierIO.write_nnkp`](@ref).
"""
function write_nnkp(filename::AbstractString, kstencil::KgridStencil; kwargs...)
    return WannierIO.write_nnkp(
        filename;
        lattice=real_lattice(reciprocal_lattice(kstencil)),
        recip_lattice=reciprocal_lattice(kstencil),
        kstencil.kgrid.kpoints,
        kstencil.kpb_k,
        kstencil.kpb_G,
        kwargs...,
    )
end
