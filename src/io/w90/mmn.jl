import WannierIO: write_mmn

"""
    $(SIGNATURES)

Write `mmn` file.

# Arguments
- `filename`: filename of the `.mmn` file
- `overlaps`: the overlap matrices
- `kstencil`: the `KSpaceStencil` struct
"""
function write_mmn(
        filename::AbstractString,
        overlaps::AbstractArray{<:Complex, 4},
        kstencil::KSpaceStencil;
        kwargs...,
    )
    return write_mmn(filename, overlaps, kstencil.kpb_k, kstencil.kpb_G; kwargs...)
end

function write_mmn(
        filename::AbstractString,
        overlaps::AbstractArray{<:Complex, 4},
        kpb_k::AbstractMatrix,
        kpb_G::AbstractMatrix;
        kwargs...,
    )
    mmn = WannierIO.Mmn(WannierIO.default_header(), overlaps, kpb_k, kpb_G)
    return WannierIO.write_mmn(filename, mmn; kwargs...)
end
