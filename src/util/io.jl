
"""
    isbinary(chars::Vector{UInt8})

Check if a sequence of chars is binary.
"""
function isbinary(chars::AbstractVector{UInt8})::Bool
    # normal ASCII chars
    text_chars = Vector{UInt8}([7, 8, 9, 10, 12, 13, 27])
    append!(text_chars, 0x20:0x99)
    deleteat!(text_chars, text_chars .== 0x7F)

    # remove normal ASCII
    filter!(x -> x ∉ text_chars, chars)

    # display([Char(_) for _ in chars])

    return length(chars) > 0
end

"""
    isbinary(filename::AbstractString)

Check if the file is in binary format.
"""
function isbinary(filename::AbstractString)
    raw_data = zeros(UInt8, 1024)

    io = open(filename)
    readbytes!(io, raw_data)
    close(io)

    return isbinary(raw_data)
end

"""
    parse_float(s::AbstractString)

Parse a string as `Float64`.

The is capable of parsing Fortran outputs, e.g. `1.0D-10`, to the ordinary `1e-10`.
"""
parse_float(s::AbstractString) = parse(Float64, replace(lowercase(strip(s)), "d" => "e"))

"""
    parse_bool(s::AbstractString)

Parse a string as `bool`.

This is capable of parsing Fortran outputs, e.g., `.true.`, `.false.`, `true`, `T`.
"""
function parse_bool(s::AbstractString)
    s = replace(lowercase(strip(s)), "." => "")[1]  # only 1st char
    return s == 't' || s == '1'
end
