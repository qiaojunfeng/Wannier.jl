
export read_spn, write_spn

"""
    read_spn(filename::AbstractString)

Read the `spn` file.

Return a `n_bands * n_bands * n_kpts * 3` array.
"""
function read_spn(filename::AbstractString)
    @info "Reading $filename"

    io = open("$filename")

    header = readline(io)

    arr = split(readline(io))
    n_bands, n_kpts = parse.(Int64, arr[1:2])

    S = zeros(ComplexF64, n_bands, n_bands, n_kpts, 3)

    for ik in 1:n_kpts
        for m in 1:n_bands
            for n in 1:m
                for s in 1:3
                    line = split(readline(io))
                    x, y = parse.(Float64, line)
                    S[n, m, ik, s] = x + im * y
                    # although each diagonal element of `S` should be real,
                    # actually it has a very small imaginary part,
                    # so we skip the conjugation on the diagonal elements.
                    m == n && continue
                    S[m, n, ik, s] = x - im * y
                end
            end
        end
    end
    @assert eof(io)
    close(io)

    println("  header  = ", header)
    println("  n_bands = ", n_bands)
    println("  n_kpts  = ", n_kpts)
    println()

    return S
end

"""
    write_spn(filename::String, S::Array{Complex,4}, header::String)

Write the `spn` file.
"""
function write_spn(
    filename::AbstractString, S::AbstractArray{T,4}, header::AbstractString
) where {T<:Complex}
    n_bands, _, n_kpts, n_spin = size(S)
    size(S, 2) == n_bands || error("S must be a square matrix")
    n_spin == 3 || error("n_spin must be 3")

    io = open(filename, "w")

    header = strip(header)
    write(io, header, "\n")

    @printf(io, "%3d %4d\n", n_bands, n_kpts)

    for ik in 1:n_kpts
        for m in 1:n_bands
            for n in 1:m
                for s in 1:3
                    a = S[n, m, ik, s]
                    @printf(io, "%26.16e  %26.16e\n", real(a), imag(a))
                end
            end
        end
    end
    close(io)

    @info "Written to file: $(filename)"
    println()

    return nothing
end

"""
    write_spn(filename::String, S::Array{Complex,4})

Write the `spn` file.

The header is "Created by Wannier.jl CURRENT_DATE".
"""
function write_spn(filename::AbstractString, S::AbstractArray{T,4}) where {T<:Complex}
    header = @sprintf "Created by Wannier.jl %s" string(now())
    return write_spn(filename, S, header)
end
