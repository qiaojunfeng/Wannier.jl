using Brillouin
using Bravais: ReciprocalBasis

export read_w90_band, write_w90_band

"""
    KPathInterpolant(kpoints, symm_idx, symm_label, recip_lattice)

Generate a `KPathInterpolant` from `kpoints` in `seedname_band.dat/kpt/labelinfo`.

# Arguments
- kpoints: fractional coordinate, each column is a kpoint.
"""
function KPathInterpolant(
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

"""
    read_w90_band(seedname::AbstractString, recip_lattice::AbstractMatrix)

# Arguments
- `recip_lattice`: each column is a reciprocal lattice vector in Cartesian coordinates.
    If given, return a tuple of `(KPathInterpolant, E)`.
    This is a more user-friendly version of
    [`read_w90_band(seedname::AbstractString)`](@ref read_w90_band(seedname::AbstractString)).

See also [`read_w90_band(seedname::AbstractString)`](@ref read_w90_band(seedname::AbstractString)).
"""
function read_w90_band(seedname::AbstractString, recip_lattice::AbstractMatrix)
    band = WannierIO.read_w90_band(seedname)
    kpi = KPathInterpolant(band.kpoints, band.symm_idx, band.symm_label, recip_lattice)
    return kpi, band.E
end

"""
    get_symm_idx_label(kpi::KPathInterpolant)

Return the symmetry indexes and labels.
"""
function get_symm_idx_label(kpi::KPathInterpolant)
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

"""
    write_w90_band(seedname::AbstractString, kpi::KPathInterpolant, E::AbstractMatrix)

Write `SEEDNAME_band.dat, SEEDNAME_band.kpt, SEEDNAME_band.labelinfo.dat`.

This is a more user-friendly version.

See also [`write_w90_band(seedname, kpoints, E, x, symm_idx, symm_label)`]
(@ref write_w90_band(seedname, kpoints, E, x, symm_idx, symm_label)).
"""
function write_w90_band(seedname::AbstractString, kpi::KPathInterpolant, E::AbstractMatrix)
    kpoints = get_kpoints(kpi::KPathInterpolant)
    x = get_x(kpi)
    symm_idx, symm_label = get_symm_idx_label(kpi)
    return WannierIO.write_w90_band(seedname, kpoints, E, x, symm_idx, symm_label)
end
