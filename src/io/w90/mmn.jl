import WannierIO: write_mmn

"""
    $(SIGNATURES)

Write `mmn` file.

# Arguments
- `filename`: filename of the `.mmn` file
- `overlaps`: the overlap matrices
- `kstencil`: the `KspaceStencil` struct
"""
function write_mmn(
    filename::AbstractString, overlaps::AbstractVector, kstencil::KspaceStencil; kwargs...
)
    return WannierIO.write_mmn(
        filename, overlaps, kstencil.kpb_k, kstencil.kpb_G; kwargs...
    )
end
