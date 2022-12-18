export read_nnkp, write_nnkp

"""
    read_nnkp(filename::AbstractString)

Read the `nnkp` file.

!!! note

    This is a wrapper of `WannierIO.read_nnkp`. It returns a [`BVectors`](@ref)
    instead of `NamedTuple`.
"""
function read_nnkp(filename::AbstractString)
    nnkp = WannierIO.read_nnkp(filename)
    n_bvecs = size(nnkp.kpb_k, 1)

    # Generate bvectors from 1st kpoint, in Cartesian coordinates
    bvectors = zeros(Float64, 3, n_bvecs)
    ik = 1
    for ib in 1:n_bvecs
        ik2 = nnkp.kpb_k[ib, ik]
        b = nnkp.kpb_b[:, ib, ik]
        bvec = nnkp.kpoints[:, ik2] + b - nnkp.kpoints[:, ik]
        bvectors[:, ib] = nnkp.recip_lattice * bvec
    end

    weights = zeros(Float64, n_bvecs)
    fill!(weights, NaN)

    return BVectors(
        nnkp.recip_lattice, nnkp.kpoints, bvectors, weights, nnkp.kpb_k, nnkp.kpb_b
    )
end

"""
    write_nnkp(filename::AbstractString, bvectors::BVectors, n_wann::Integer)

Write nnkp that can be used by `pw2wannier90`.

!!! note

    This is a wrapper of `WannierIO.write_nnkp`.
"""
function write_nnkp(
    filename::AbstractString,
    bvectors::BVectors,
    n_wann::Integer,
    exclude_bands::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)
    return WannierIO.write_nnkp(
        filename,
        bvectors.recip_lattice,
        bvectors.kpoints,
        bvectors.kpb_k,
        bvectors.kpb_b,
        n_wann,
        exclude_bands,
    )
end
