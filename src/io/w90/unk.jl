using Printf: @printf, @sprintf

@doc raw"""
Read UNK file for Bloch wavefunctions.
"""
function read_unk(filename::String)
    @info "Reading unk file:" filename

    io = open(filename)

    line = split(strip(readline(io)))
    n_gx, n_gy, n_gz, ik, n_bands = parse.(Int, line)

    Ψ = zeros(ComplexF64, n_gx, n_gy, n_gz, n_bands)

    for ib in 1:n_bands
        for iz in 1:n_gz
            for iy in 1:n_gy
                for ix in 1:n_gx
                    line = split(strip(readline(io)))
                    v1, v2 = parse.(Float64, line)
                    Ψ[ix, iy, iz, ib] = v1 + im * v2
                end
            end
        end
    end

    close(io)

    println("  n_gx    = ", n_gx)
    println("  n_gy    = ", n_gy)
    println("  n_gz    = ", n_gz)
    println("  ik      = ", ik)
    println("  n_bands = ", n_bands)
    println()

    # ik: at which kpoint? start from 1
    return ik, Ψ
end

@doc raw"""
Write UNK file for Bloch wavefunctions.
    ik: at which kpoint? start from 1
"""
function write_unk(filename::String, ik::Int, Ψ::Array{T,4}) where {T<:Complex}
    @info "Writing unk file:" filename

    n_gx, n_gy, n_gz, n_bands = size(Ψ)

    open(filename, "w") do io
        @printf(io, " %11d %11d %11d %11d %11d\n", n_gx, n_gy, n_gz, ik, n_bands)

        for ib in 1:n_bands
            for iz in 1:n_gz
                for iy in 1:n_gy
                    for ix in 1:n_gx
                        v1 = real(Ψ[ix, iy, iz, ib])
                        v2 = imag(Ψ[ix, iy, iz, ib])
                        @printf(io, " %18.10e %18.10e\n", v1, v2)
                    end
                end
            end
        end
    end

    @info "Written to file: $(filename)"
    println("  n_gx    = ", n_gx)
    println("  n_gy    = ", n_gy)
    println("  n_gz    = ", n_gz)
    println("  ik      = ", ik)
    println("  n_bands = ", n_bands)
    println()

    return nothing
end
