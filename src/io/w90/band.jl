using Printf: @printf, @sprintf
import DelimitedFiles as Dlm

"""
read si_band.dat  si_band.kpt  si_band.labelinfo.dat
"""
function read_w90_band(seedname::String)
    kpoints = Dlm.readdlm("$(seedname)_band.kpt", Float64; skipstart=1)
    # remove weights, then transpose, last idx is kpt
    kpoints = Matrix(transpose(kpoints[:, 1:3]))
    n_kpts = size(kpoints, 2)

    dat = Dlm.readdlm("$(seedname)_band.dat", Float64)
    x = reshape(dat[:, 1], n_kpts, :)[:, 1]
    E = reshape(dat[:, 2], n_kpts, :)
    E = Matrix(transpose(E))
    n_bands = size(E, 1)

    io = open("$(seedname)_band.labelinfo.dat")
    labels = readlines(io)
    close(io)

    n_symm = length(labels)
    symm_idx = Vector{Int}(undef, n_symm)
    symm_label = Vector{String}(undef, n_symm)
    for (i, line) in enumerate(labels)
        lab, idx = split(line)[1:2]
        symm_idx[i] = parse(Int, idx)
        symm_label[i] = lab
    end

    return (kpoints=kpoints, E=E, x=x, symm_idx=symm_idx, symm_label=symm_label)
end

"""
write si_band.dat  si_band.kpt  si_band.labelinfo.dat
"""
function write_w90_band(
    seedname::String,
    kpoints::AbstractMatrix{T},
    E::AbstractMatrix{T},
    x::AbstractVector{T},
    symm_idx::AbstractVector{Int},
    symm_label::AbstractVector{String},
) where {T<:Real}
    n_kpts = size(kpoints, 2)
    size(kpoints, 1) != 3 && error("kpoints must be 3 x n_kpts")

    n_bands = size(E, 1)
    size(E, 2) != n_kpts && error("E must be n_bands x n_kpts")

    length(x) != n_kpts && error("x must be n_kpts")

    n_symm = length(symm_idx)
    n_symm != length(symm_label) && error("symm_idx and symm_label must be same length")

    filename = "$(seedname)_band.kpt"
    open(filename, "w") do io
        @printf(io, "       %5d\n", n_kpts)
        for ik in 1:n_kpts
            @printf(io, "  %10.6f  %10.6f  %10.6f   1.0\n", kpoints[:, ik]...)
        end
    end
    @info "Written to $filename"

    filename = "$(seedname)_band.dat"
    open(filename, "w") do io
        for ib in 1:n_bands
            for ik in 1:n_kpts
                @printf(io, " %15.8E %15.8E\n", x[ik], E[ib, ik])
            end
        end
    end
    @info "Written to $filename"

    filename = "$(seedname)_band.labelinfo.dat"
    open(filename, "w") do io
        for i in 1:n_symm
            idx = symm_idx[i]
            @printf(
                io,
                "%2s %31d %20.10f %17.10f %17.10f %17.10f\n",
                symm_label[i],
                idx,
                x[i],
                kpoints[:, idx]...
            )
        end
    end
    @info "Written to $filename"
    println()
    return nothing
end
