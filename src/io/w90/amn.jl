using Printf: @printf, @sprintf
using Dates: now

function read_amn(filename::String)
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
Wrapper function to read AMN and Lowdin orthonormalize it.

Usually A matrix needs to be unitary or semi-unitary.
"""
function read_orthonorm_amn(filename::String)
    A = read_amn(filename)
    A .= orthonorm_lowdin(A)
    return A
end

"""
Output amn file
"""
function write_amn(filename::String, A::Array{ComplexF64,3}, header::String)
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

"""
Output amn file
"""
function write_amn(filename::String, A::Array{ComplexF64,3})
    header = @sprintf "Created by Wannier.jl %s" string(now())
    return write_amn(filename, A, header)
end
