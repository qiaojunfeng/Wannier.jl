using Printf: @printf, @sprintf
using Dates: now

export read_amn, read_orthonorm_amn, write_amn

"""
    read_amn(filename::AbstractString)

Read the `amn` file.

Return a `n_bands * n_wann * n_kpts` array.

See also [`read_orthonorm_amn`](@ref read_orthonorm_amn).
"""
function read_amn(filename::AbstractString)
    @info "Reading $filename"

    io = open("$filename")

    header = readline(io)

    arr = split(readline(io))
    n_bands, n_kpts, n_wann = parse.(Int64, arr[1:3])

    A = zeros(ComplexF64, n_bands, n_wann, n_kpts)

    while !eof(io)
        line = readline(io)
        arr = split(line)

        m, n, k = parse.(Int64, arr[1:3])

        a = parse(Float64, arr[4]) + im * parse(Float64, arr[5])
        A[m, n, k] = a
    end

    close(io)

    println("  header  = ", header)
    println("  n_bands = ", n_bands)
    println("  n_wann  = ", n_wann)
    println("  n_kpts  = ", n_kpts)
    println()

    return A
end

"""
    read_orthonorm_amn(filename::AbstractString)

Read and orthonormalize the `amn` file.

Wrapper function to read `amn` and Lowdin orthonormalize it.
The `A` matrix needs to be unitary or semi-unitary,
so in most cases this function should be used instead of [`read_amn`](@ref read_amn).

See also [`read_amn`](@ref read_amn).
"""
function read_orthonorm_amn(filename::AbstractString)
    A = read_amn(filename)
    A .= orthonorm_lowdin(A)
    return A
end

"""
    write_amn(filename::AbstractString, A::Array{ComplexF64,3})
    write_amn(filename::AbstractString, A::Array{ComplexF64,3}, header::AbstractString)

Output `amn` file.

# Arguments
- `header`: optional, default is "Created by Wannier.jl CURRENT_DATE"
"""
function write_amn(filename::AbstractString, A::Array{ComplexF64,3}, header::AbstractString)
    n_bands, n_wann, n_kpts = size(A)

    io = open(filename, "w")

    header = strip(header)
    write(io, header, "\n")

    @printf(io, "%3d %4d %4d\n", n_bands, n_kpts, n_wann)

    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_bands
                a = A[ib, iw, ik]
                @printf(io, "%5d %4d %4d  %16.12f  %16.12f\n", ib, iw, ik, real(a), imag(a))
            end
        end
    end

    close(io)

    @info "Written to file: $(filename)"
    println()

    return nothing
end

function write_amn(filename::String, A::Array{ComplexF64,3})
    header = @sprintf "Created by Wannier.jl %s" string(now())
    return write_amn(filename, A, header)
end
