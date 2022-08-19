using Printf: @printf, @sprintf
import DelimitedFiles as Dlm
using Brillouin
using Bravais: ReciprocalBasis

function _read_w90_band_kpt(filename::AbstractString)
    # in fractional coordinates
    kpoints = Dlm.readdlm(filename, Float64; skipstart=1)
    weights = kpoints[:, 4]
    # remove weights, then transpose, last idx is kpt
    kpoints = Matrix(transpose(kpoints[:, 1:3]))
    return kpoints, weights
end

function _read_w90_band_dat(filename::AbstractString)
    # 1st read to get n_kpts
    io = open(filename)
    n_kpts = 0
    while true
        line = strip(readline(io))
        isempty(line) && break
        n_kpts += 1
    end

    seekstart(io)
    dat = Dlm.readdlm(io, Float64)
    close(io)

    x = reshape(dat[:, 1], n_kpts, :)[:, 1]
    E = reshape(dat[:, 2], n_kpts, :)
    # first dim is n_bands, 2nd dim is n_kpts
    E = Matrix(transpose(E))
    return x, E
end

function _read_w90_band_labelinfo(filename::AbstractString)
    io = open(filename)
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

    return symm_idx, symm_label
end

"""
read si_band.dat  si_band.kpt  si_band.labelinfo.dat
"""
function read_w90_band(seedname::AbstractString)
    kpoints, _ = _read_w90_band_kpt("$(seedname)_band.kpt")
    x, E = _read_w90_band_dat("$(seedname)_band.dat")
    symm_idx, symm_label = _read_w90_band_labelinfo("$(seedname)_band.labelinfo.dat")
    return (kpoints=kpoints, E=E, x=x, symm_idx=symm_idx, symm_label=symm_label)
end

"""
Generate a KPathInterpolant from kpoints in seedname_band.dat/kpt/labelinfo.

kpoints: fractional coordinate
"""
function _get_kpath_interpolant(
    kpoints::AbstractMatrix,
    symm_idx::AbstractVector{T},
    symm_label::AbstractVector{R},
    recip_lattice::AbstractMatrix,
) where {T<:Integer,R<:AbstractString}
    # kpoints along path
    kpaths = Vector{Vector{Vec3{Float64}}}()
    # symmetry points
    labels = Vector{Dict{Int,Symbol}}()

    i0 = symm_idx[1]  # 1st point
    lab = Dict{Int,Symbol}()  # label of each line
    push!(lab, i0 => Symbol(symm_label[1]))
    for (i, l) in zip(symm_idx[2:end], symm_label[2:end])
        if i == i0 + 1
            push!(labels, lab)
            lab = Dict{Int,Symbol}()
        end
        push!(lab, i => Symbol(l))
        i0 = i
    end
    push!(labels, lab)

    for lab in labels
        kp = Vector{Vec3{Float64}}()  # kpath of each line
        ik1 = minimum(keys(lab))
        ik2 = maximum(keys(lab))
        append!(kp, [v for v in eachcol(kpoints[:, ik1:ik2])])
        push!(kpaths, kp)
    end

    for (i, lab) in enumerate(labels)
        ik1 = minimum(keys(lab)) - 1
        labels[i] = Dict(((k - ik1) => v) for (k, v) in lab)
    end

    basis = ReciprocalBasis([v for v in eachcol(recip_lattice)])
    setting = Ref(Brillouin.LATTICE)
    kpi = KPathInterpolant(kpaths, labels, basis, setting)
    return kpi
end

function read_w90_band(seedname::AbstractString, recip_lattice::AbstractMatrix)
    band = read_w90_band(seedname)
    kpi = _get_kpath_interpolant(
        band.kpoints, band.symm_idx, band.symm_label, recip_lattice
    )
    return kpi, band.E
end

function _write_w90_band_kpt(
    filename::AbstractString,
    kpoints::AbstractMatrix{T},
    weights::Union{Nothing,AbstractVector{T}}=nothing,
) where {T<:Real}
    n_kpts = size(kpoints, 2)
    size(kpoints, 1) == 3 || error("kpoints must be 3 x n_kpts")

    isnothing(weights) && (weights = ones(T, n_kpts))
    length(weights) == n_kpts || error("weights must be n_kpts")

    open(filename, "w") do io
        @printf(io, "       %5d\n", n_kpts)
        for ik in 1:n_kpts
            k = kpoints[:, ik]
            w = weights[ik]
            @printf(io, "  %10.6f  %10.6f  %10.6f   %10.6f\n", k..., w)
        end
    end
    @info "Written to $filename"
end

function _write_w90_band_dat(
    filename::AbstractString, x::AbstractVector{T}, E::AbstractMatrix{T}
) where {T<:Real}
    n_bands, n_kpts = size(E)
    length(x) == n_kpts || error("x must be n_kpts")

    open(filename, "w") do io
        for ib in 1:n_bands
            for ik in 1:n_kpts
                @printf(io, " %15.8E %15.8E\n", x[ik], E[ib, ik])
            end
            @printf(io, "\n")
        end
    end
    @info "Written to $filename"
end

function _write_w90_band_labelinfo(
    filename::AbstractString,
    symm_idx::AbstractVector{T},
    symm_label::AbstractVector{R},
    x::AbstractVector{S},
    kpoints::AbstractMatrix{S},
) where {T<:Integer,R<:AbstractString,S<:Real}
    n_symm = length(symm_idx)
    n_symm == length(symm_label) || error("symm_idx and symm_label must be same length")

    open(filename, "w") do io
        for i in 1:n_symm
            idx = symm_idx[i]
            @printf(
                io,
                "%2s %31d %20.10f %17.10f %17.10f %17.10f\n",
                symm_label[i],
                idx,
                x[idx],
                kpoints[:, idx]...
            )
        end
    end
    @info "Written to $filename"
end

"""
write si_band.dat  si_band.kpt  si_band.labelinfo.dat
"""
function write_w90_band(
    seedname::AbstractString,
    kpoints::AbstractMatrix{T},
    E::AbstractMatrix{T},
    x::AbstractVector{T},
    symm_idx::AbstractVector{R},
    symm_label::AbstractVector{S},
) where {T<:Real,R<:Integer,S<:AbstractString}
    size(kpoints, 2) == size(E, 2) || error("kpoints and E have different n_kpts")

    filename = "$(seedname)_band.kpt"
    _write_w90_band_kpt(filename, kpoints)

    filename = "$(seedname)_band.dat"
    _write_w90_band_dat(filename, x, E)

    filename = "$(seedname)_band.labelinfo.dat"
    _write_w90_band_labelinfo(filename, symm_idx, symm_label, x, kpoints)

    println()
    return nothing
end

function get_kpoints(kpi::KPathInterpolant)
    kpi_frac = latticize(kpi)

    n_kpts = sum(length(kp) for kp in kpi_frac.kpaths)
    kpoints = zeros(Float64, 3, n_kpts)

    ik = 1
    for kp in kpi_frac.kpaths
        for k in kp
            kpoints[:, ik] = k
            ik += 1
        end
    end
    return kpoints
end

function _get_symm_idx_label(kpi::KPathInterpolant)
    kpi_frac = latticize(kpi)

    symm_idx = Vector{Int}()
    symm_label = Vector{String}()
    ik0 = 0
    for lab in kpi_frac.labels
        for (ik, l) in lab
            push!(symm_idx, ik + ik0)
            push!(symm_label, String(l))
        end
        ik0 = maximum(keys(lab))
    end

    perm = sortperm(symm_idx)
    symm_idx = symm_idx[perm]
    symm_label = symm_label[perm]

    return symm_idx, symm_label
end

function write_w90_band(seedname::AbstractString, kpi::KPathInterpolant, E::AbstractMatrix)
    kpoints = get_kpoints(kpi::KPathInterpolant)
    x = get_x(kpi)
    symm_idx, symm_label = _get_symm_idx_label(kpi)
    return write_w90_band(seedname, kpoints, E, x, symm_idx, symm_label)
end
