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
    return KspaceStencil(nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G)
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

"""
    $(SIGNATURES)

Check if the kpoint stencil (bvectors) contain 6 cubic neighbors.

[`parallel_transport`](@ref) requires the 6 nearest neighbors.

# Arguments
- `kpoints`: the kpoints in fractional coordinates
- `kpb_k`: indices of kpoint (in `kpoints`) that is translational equivalent to k+b
- `kpb_G`: translation vector that maps k to k+b
"""
function has_cubic_neighbors(
    kpoints::AbstractVector,
    kpb_k::AbstractVector,
    kpb_G::AbstractVector;
    kgrid_size::AbstractVector=guess_kgrid_size(kpoints),
    atol::AbstractFloat=1e-6,
)
    dkx, dky, dkz = 1 ./ kgrid_size
    # In fractional coordinates
    bvectors_cubic = [
        [dkx, 0.0, 0.0],
        [0.0, dky, 0.0],
        [0.0, 0.0, dkz],
        [-dkx, 0.0, 0.0],
        [0.0, -dky, 0.0],
        [0.0, 0.0, -dkz],
    ]

    # Only need to check the first kpoint
    idx = 1
    bvectors = kpoints[kpb_k[idx]] .+ kpb_G[idx] .- Ref(kpoints[idx])

    for b in bvectors_cubic
        if isnothing(findfirst(isapprox(b; atol), bvectors))
            return false
        end
    end
    return true
end

function has_cubic_neighbors(kstencil::KspaceStencil; atol::AbstractFloat=1e-6)
    return has_cubic_neighbors(
        kstencil.kpoints, kstencil.kpb_k, kstencil.kpb_G; kstencil.kgrid_size, atol
    )
end

function has_cubic_neighbors(filename::AbstractString; atol::AbstractFloat=1e-6)
    nnkp = read_nnkp(filename)
    return has_cubic_neighbors(nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G; atol)
end
