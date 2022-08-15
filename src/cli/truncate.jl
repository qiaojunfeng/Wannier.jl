"""
Truncate MMN/EIG/UNK files.

# Args

- `seedname`: seedname for MMN/EIG files
- `keep-bands`: indices of bands to be kept, start from 1

# Options

- `--outdir`: dirname for output MMN/EIG. Default is `truncate`

# Flags

- `--rotate-unk`: also truncate UNK files, for plotting WFs
"""
@cast function trunc(
    seedname::String, keep_bands::Int...; outdir::String="truncate", rotate_unk::Bool=false
)
    # tuple to vector
    keep_bands = collect(keep_bands)
    truncate_w90(seedname, keep_bands, outdir, rotate_unk)
    return nothing
end
