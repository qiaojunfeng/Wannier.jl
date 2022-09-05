"""
Truncate number of bands in `mmn`/`eig`/`unk` files.

# Args

- `seedname`: seedname for `mmn`/`eig` files
- `keep-bands`: indices of bands to be kept, start from 1

# Options

- `--outdir`: dirname for output `mmn`/`eig`. Default is `truncate`

# Flags

- `--rotate-unk`: also truncate `unk` files, for plotting WFs
"""
@cast function trunc(
    seedname::String, keep_bands::Int...; outdir::String="truncate", rotate_unk::Bool=false
)
    # tuple to vector
    keep_bands = collect(keep_bands)
    truncate_w90(seedname, keep_bands, outdir, rotate_unk)
    return nothing
end
