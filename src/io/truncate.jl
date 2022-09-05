export truncate_w90

"""
    truncate_mmn_eig(seedname, keep_bands::Vector{Int}, outdir="truncate")

Truncate number of bands of `mmn` and `eig` files.

# Arguments
- `keep_bands`: a vector of band indices to keep, starting from 1.
- `outdir`: the folder for writing `mmn` and `eig` files.

!!! tip

    This is useful for generating valence only `mmn`, `eig` files from a
    valence+conduction NSCF calculation, so that no need to recompute NSCF with
    lower number of bands again.
"""
function truncate_mmn_eig(
    seedname::AbstractString,
    keep_bands::AbstractVector{Int},
    outdir::AbstractString="truncate",
)
    !isdir(outdir) && mkdir(outdir)

    # for safety, in case seedname = "../si" then joinpath(outdir, seedname)
    # will overwrite the original file
    seedname_base = basename(seedname)

    E = read_eig("$seedname.eig")
    E1 = E[keep_bands, :]
    write_eig(joinpath(outdir, "$seedname_base.eig"), E1)

    M, kpb_k, kpb_b = read_mmn("$seedname.mmn")
    M1 = M[keep_bands, keep_bands, :, :]
    write_mmn(joinpath(outdir, "$seedname_base.mmn"), M1, kpb_k, kpb_b)

    return nothing
end

"""
    truncate_unk(dir, keep_bands::Vector{Int}, outdir="truncate")

Truncate `UNK` files for specified bands.

# Arguments
- `dir`: folder of `UNK` files.
- `keep_bands`: the band indexes to keep. Start from 1.
- `outdir`: folder to write output `UNK` files.
"""
function truncate_unk(
    dir::AbstractString, keep_bands::AbstractVector{Int}, outdir::AbstractString="truncate"
)
    !isdir(outdir) && mkdir(outdir)

    regex = r"UNK(\d{5})\.\d"

    for unk in readdir(dir)
        m = match(regex, unk)
        m === nothing && continue

        println(unk)
        # for safety, in case unk = "../UNK00001.1" then joinpath(dir, unk)
        # will overwrite the original file
        unk_base = basename(unk)

        ik = parse(Int, m.captures[1])
        ik1, Ψ = read_unk(joinpath(dir, unk))
        @assert ik == ik1

        Ψ1 = Ψ[:, :, :, keep_bands]
        write_unk(joinpath(outdir, unk), ik, Ψ1)
    end

    return nothing
end

"""
    truncate_w90(seedname, keep_bands::Vector{Int}, outdir="truncate", unk=false)

Truncate `mmn`, `eig`, and optionally `UNK` files.

# Arguments
- seedname: seedname for input `mmn` and `eig` files.
- keep_bands: Band indexes to be kept, start from 1.
- unk: If true also truncate `UNK` files.
- outdir: folder for output files.
"""
function truncate_w90(
    seedname::AbstractString,
    keep_bands::AbstractVector{Int},
    outdir::AbstractString="truncate",
    unk::Bool=false,
)
    @info "Truncat AMN/MMN/EIG files"

    !isdir(outdir) && mkdir(outdir)

    # E = read_eig("$seedname.eig")
    # n_bands = size(E, 1)
    # keep_bands = [i for i = 1:n_bands if i ∉ exclude_bands]

    truncate_mmn_eig(seedname, keep_bands, outdir)

    dir = dirname(seedname)

    if unk
        truncate_unk(dir, keep_bands, outdir)
    end

    println("Truncated files written in ", outdir)
    return nothing
end

"""
    truncate(model::Model, keep_bands::Vector{Int}, keep_wfs::Vector{Int}=nothing;
        orthonorm_A::Bool=true)

Truncate `A`, `M`, `E` matrices in `model`.

# Arguments
- `model`: the `Model` to be truncated.
- `keep_bands`: Band indexes to be kept, start from 1.
- `keep_wfs`: WF indexes to be kept, start from 1. If `nothing`, keep all.

# Keyword arguments
- `orthonorm_A`: If true, Lowdin orthonormalize `A` after truncation.
    The `A` needs to be (semi-)unitary, so it should always be true.
"""
function truncate(
    model::Model, keep_bands::T, keep_wfs::Union{T,Nothing}=nothing; orthonorm_A::Bool=true
) where {T<:AbstractVector{Int}}
    all(1 .<= keep_bands .<= model.n_bands) || error("Invalid band index")
    if !isnothing(keep_wfs)
        all(1 .<= keep_wfs .<= model.n_wann) || error("Invalid WF index")
        length(keep_wfs) <= length(keep_bands) || error("Number of WFs > number of bands")
    end

    E = model.E[keep_bands, :]
    M = model.M[keep_bands, keep_bands, :, :]
    A = model.A[keep_bands, :, :]

    if !isnothing(keep_wfs)
        A = A[:, keep_wfs, :]
    end
    if orthonorm_A
        A = orthonorm_lowdin(A)
    end
    frozen_bands = model.frozen_bands[keep_bands, :]

    model2 = Model(
        model.lattice,
        model.atom_positions,
        model.atom_labels,
        model.kgrid,
        model.kpoints,
        model.bvectors,
        frozen_bands,
        M,
        A,
        E,
    )
    return model2
end
