using Brillouin
using Bravais: ReciprocalBasis

export read_w90_band, write_w90_band

"""
    KPathInterpolant(kpoints, symm_point_indices, symm_point_labels, recip_lattice)

Generate a `KPathInterpolant` from `kpoints` in `seedname_band.dat/kpt/labelinfo`.

# Arguments
- kpoints: fractional coordinate, each column is a kpoint.
"""
function KPathInterpolant(
    kpoints::AbstractVector,
    symm_point_indices::AbstractVector{T},
    symm_point_labels::AbstractVector{R},
    recip_lattice::AbstractMatrix,
) where {T<:Integer,R<:AbstractString}
    # kpoints along path
    kpaths = Vector{Vector{Vec3{Float64}}}()
    # symmetry points
    labels = Vector{Dict{Int,Symbol}}()

    i0 = symm_point_indices[1]  # 1st point
    lab = Dict{Int,Symbol}()  # label of each line
    push!(lab, i0 => Symbol(symm_point_labels[1]))
    for (i, l) in zip(symm_point_indices[2:end], symm_point_labels[2:end])
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
        append!(kp, [v for v in kpoints[ik1:ik2]])
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
    kpi = KPathInterpolant(
        band.kpoints, band.symm_point_indices, band.symm_point_labels, recip_lattice
    )
    return kpi, band.eigenvalues
end

"""
    get_symm_point_indices_label(kpi::KPathInterpolant)

Return the symmetry indexes and labels.
"""
function get_symm_point_indices_label(kpi::KPathInterpolant)
    kpi_frac = latticize(kpi)

    symm_point_indices = Vector{Int}()
    symm_point_labels = Vector{String}()
    ik0 = 0
    for lab in kpi_frac.labels
        for (ik, l) in lab
            push!(symm_point_indices, ik + ik0)
            push!(symm_point_labels, String(l))
        end
        ik0 += maximum(keys(lab))
    end

    perm = sortperm(symm_point_indices)
    symm_point_indices = symm_point_indices[perm]
    symm_point_labels = symm_point_labels[perm]

    return symm_point_indices, symm_point_labels
end

"""
    write_w90_band(seedname::AbstractString, kpi::KPathInterpolant, E::AbstractMatrix)

Write `SEEDNAME_band.dat, SEEDNAME_band.kpt, SEEDNAME_band.labelinfo.dat`.

This is a more user-friendly version.

See also [`write_w90_band(seedname, kpoints, E, x, symm_point_indices, symm_point_labels)`]
(@ref write_w90_band(seedname, kpoints, E, x, symm_point_indices, symm_point_labels)).
"""
function write_w90_band(prefix::AbstractString, kpi::KPathInterpolant, eigenvalues::Vector)
    kpoints = get_kpoints(kpi::KPathInterpolant)
    x = get_x(kpi)
    symm_point_indices, symm_point_labels = get_symm_point_indices_label(kpi)
    return WannierIO.write_w90_band(
        prefix; x, eigenvalues, kpoints, symm_point_indices, symm_point_labels
    )
end

"""
    write_w90_kpt_label(seedname::AbstractString, kpi::KPathInterpolant)

Write `SEEDNAME_band.kpt` and `SEEDNAME_band.labelinfo.dat`.

This allows generating the high-symmetry kpoints and labels from crystal
structure, and use the generated kpoints in `pw.x` `bands` calculation
or in the `win` input file for `Wannier90`.

# Example
```julia
win = read_win("si2.win")
kp = get_kpath(win.unit_cell_cart, win.atoms_frac, win.atom_labels)
kpi = Wannier.interpolate_w90(kp, 100)
Wannier.write_w90_kpt_label("si2", kpi)
```
"""
function write_w90_kpt_label(prefix::AbstractString, kpi::KPathInterpolant)
    kpoints = get_kpoints(kpi::KPathInterpolant)
    x = get_x(kpi)
    symm_point_indices, symm_point_labels = get_symm_point_indices_label(kpi)

    filename = "$(prefix)_band.kpt"
    WannierIO.write_w90_band_kpt(filename; kpoints)

    filename = "$(prefix)_band.labelinfo.dat"
    WannierIO.write_w90_band_labelinfo(
        filename; x, kpoints, symm_point_indices, symm_point_labels
    )

    println()
    return nothing
end
