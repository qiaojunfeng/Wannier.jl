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
    return KspaceStencil(nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"])
end

"""
    $(SIGNATURES)

# Arguments
- `filename`: the filename of the `nnkp` file
- `weights`: the weights of the b-vectors. If `nothing`, the weights are set to 0.0.
    This is useful for [`parallel_transport`](@ref) using 6 cubic neighbors, where
    the the ``b``-vectors are not complete, thus not possible to compute the weights.
"""
function WannierIO.read_nnkp(
        filename::AbstractString, weights::Union{AbstractVector, Nothing}
    )
    nnkp = WannierIO.read_nnkp(filename)

    kpoints = nnkp["kpoints"]
    recip_lattice = nnkp["recip_lattice"]
    kpb_k = nnkp["kpb_k"]
    kpb_G = nnkp["kpb_G"]
    kgrid_size = guess_kgrid_size(kpoints)
    bvectors = get_bvectors(recip_lattice, kpoints, kpb_k, kpb_G)

    n_bvecs = length(kpb_k[1])
    if isnothing(weights)
        # If weights are not provided, set them to 0.0
        weights = zeros(Float64, n_bvecs)
    else
        (length(weights) == n_bvecs) ||
            error("length of weights does not match number of b-vectors")
    end

    return KspaceStencil{Float64}(
        recip_lattice, kgrid_size, kpoints, bvectors, weights, kpb_k, kpb_G
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
    params = OrderedDict{String, Any}(
        "lattice" => real_lattice(kstencil),
        "recip_lattice" => reciprocal_lattice(kstencil),
        "kpoints" => kstencil.kpoints,
        "kpb_k" => kstencil.kpb_k,
        "kpb_G" => kstencil.kpb_G,
    )
    for (k, v) in kwargs
        params[string(k)] = v
    end
    return WannierIO.write_nnkp(filename, params)
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
        kgrid_size::AbstractVector = guess_kgrid_size(kpoints),
        atol::AbstractFloat = 1.0e-6,
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

function has_cubic_neighbors(kstencil::KspaceStencil; atol::AbstractFloat = 1.0e-6)
    return has_cubic_neighbors(
        kstencil.kpoints, kstencil.kpb_k, kstencil.kpb_G; kstencil.kgrid_size, atol
    )
end

function has_cubic_neighbors(filename::AbstractString; atol::AbstractFloat = 1.0e-6)
    nnkp = read_nnkp(filename)
    return has_cubic_neighbors(nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]; atol)
end

"""
    $(SIGNATURES)

Write a nnkp file with 6 cubic neighbors. Useful for [`parallel_transport`](@ref).

!!! note

    `nnkp` file also contains the `exclude_bands` parameters, which will affect
    the `mmn` file computed by `pw2wannier90.x`.

# Arguments
- `filename`: the filename of the new `nnkp` file
- `win`: the Wannier90 parameters, e.g. returned by [`read_win`](@ref)
"""
function write_nnkp_cubic(filename::AbstractString, win::Union{NamedTuple, AbstractDict})
    recip_latt = reciprocal_lattice(win["unit_cell_cart"])
    kstencil = generate_kspace_stencil(
        recip_latt, win["mp_grid"], win["kpoints"], Wannier.CubicNearestKspaceStencil()
    )
    return write_nnkp(
        filename,
        kstencil;
        exclude_bands = get(win, "exclude_bands", nothing),
        # Need a fake projections block such that pw2wannier90.x can run
        projections = WannierIO.HydrogenOrbital[],
    )
end

function write_nnkp_cubic(filename::AbstractString, win::AbstractString)
    return write_nnkp_cubic(filename, read_win(win))
end
