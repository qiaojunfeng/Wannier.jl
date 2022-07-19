using Printf: @printf, @sprintf
using Dates: now

function read_mmn(filename::String)
    @info "Reading mmn file: $filename"
    io = open(filename)

    # skip header
    header = readline(io)

    line = readline(io)
    n_bands, n_kpts, n_bvecs = map(x -> parse(Int64, x), split(line))

    # overlap matrix
    M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    # for each point, list of neighbors, (K) representation
    kpb_k = zeros(Int64, n_bvecs, n_kpts)
    kpb_b = zeros(Int64, 3, n_bvecs, n_kpts)

    while !eof(io)
        for ib in 1:n_bvecs
            line = readline(io)
            arr = split(line)
            k = parse(Int64, arr[1])
            kpb_k[ib, k] = parse(Int64, arr[2])
            kpb_b[:, ib, k] = map(x -> parse(Int64, x), arr[3:5])
            for n in 1:n_bands
                for m in 1:n_bands
                    line = readline(io)
                    arr = split(line)
                    o = parse(Float64, arr[1]) + im * parse(Float64, arr[2])
                    M[m, n, ib, k] = o
                end
            end
        end
    end

    @assert !any(isnan.(M))

    close(io)

    println("  header  = ", header)
    println("  n_bands = ", n_bands)
    println("  n_bvecs = ", n_bvecs)
    println("  n_kpts  = ", n_kpts)
    println()

    return M, kpb_k, kpb_b
end

"""
Output mmn file
"""
function write_mmn(
    filename::String,
    M::Array{ComplexF64,4},
    kpb_k::Matrix{Int},
    kpb_b::Array{Int,3},
    header::String,
)
    n_bands, _, n_bvecs, n_kpts = size(M)
    n_bands != size(M, 2) && error("M must be n_bands x n_bands x n_bvecs x n_kpts")

    size(kpb_k) != (n_bvecs, n_kpts) && error("kpb_k has wrong size")
    size(kpb_b) != (3, n_bvecs, n_kpts) && error("kpb_b has wrong size")

    open(filename, "w") do io
        header = strip(header)
        write(io, header, "\n")

        @printf(io, "    %d   %d    %d \n", n_bands, n_kpts, n_bvecs)

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                @printf(io, "%d %d %d %d %d\n", ik, kpb_k[ib, ik], kpb_b[:, ib, ik]...)

                for n in 1:n_bands
                    for m in 1:n_bands
                        o = M[m, n, ib, ik]
                        @printf(io, "  %16.12f  %16.12f\n", real(o), imag(o))
                    end
                end
            end
        end
    end

    @info "Written to file: $(filename)"
    println()

    return nothing
end

"""
Output mmn file
"""
function write_mmn(
    filename::String, M::Array{ComplexF64,4}, kpb_k::Matrix{Int}, kpb_b::Array{Int,3}
)
    header = @sprintf "Created by Wannier.jl %s" string(now())
    return write_mmn(filename, M, kpb_k, kpb_b, header)
end
