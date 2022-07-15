
"""Check if string is binary."""
function isbinary(chars::Vector{UInt8})::Bool
    # normal ASCII chars
    text_chars = Vector{UInt8}([7, 8, 9, 10, 12, 13, 27])
    append!(text_chars, 0x20:0x99)
    deleteat!(text_chars, text_chars .== 0x7F)

    # remove normal ASCII
    filter!(x -> x ∉ text_chars, chars)

    # display([Char(_) for _ in chars])

    return length(chars) > 0
end

"""Check if file is in binary format"""
function isbinary_file(filename::String)
    raw_data = zeros(UInt8, 1024)

    io = open(filename)
    readbytes!(io, raw_data)
    close(io)

    return isbinary(raw_data)
end

"""
Parse a string as Float64.

Fortran some times use e.g. 1.0D-10 for 1e-10.
"""
parse_float(str::AbstractString) = parse(Float64, replace(lowercase(str), "d" => "e"))
