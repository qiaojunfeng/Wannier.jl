using Bravais: ReciprocalBasis
using Brillouin: LATTICE, KPathInterpolant
# to extend some methods
import WannierIO: read_w90_band, write_w90_band

export read_w90_band, write_w90_band

"""
    $(SIGNATURES)

Generate a `KPathInterpolant` from the kpoint coordinates and high-symmetry kpoints.

The `WannierIO.read_w90_band(prefix)` function returns the required arguments
of this function.

# Arguments
- `recip_lattice`: each column is a reciprocal lattice vector in Å
- `kpoints`: list of kpoint coordinates along a kpath, fractional coordinates
- `symm_point_indices`: indices of the high-symmetry kpoints
- `symm_point_labels`: labels of the high-symmetry kpoints
"""
function generate_w90_kpoint_path(
    recip_lattice::AbstractMatrix,
    kpoints::AbstractVector,
    symm_point_indices::AbstractVector,
    symm_point_labels::AbstractVector,
)
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
    setting = Ref(LATTICE)
    kpi = KPathInterpolant(kpaths, labels, basis, setting)
    return kpi
end

"""
    $(SIGNATURES)

# Arguments
- `recip_lattice`: each column is a reciprocal lattice vector in Å.

# Return
- `kpi`: a `KPathInterpolant` object
- `eigenvalues`: the eigenvalues of the band structure

!!! note

    The `WannierIO.read_w90_band(prefix)` function returns a `NamedTuple`
    containing basis variables such as `kpoints`, `symm_point_indices`, etc.
    Here, once we know the `recip_lattice`, we can generate a `KPathInterpolant`
    which can be used directly in plotting functions.
"""
function read_w90_band(prefix::AbstractString, recip_lattice::AbstractMatrix)
    band = WannierIO.read_w90_band(prefix)
    kpi = generate_w90_kpoint_path(
        recip_lattice, band.kpoints, band.symm_point_indices, band.symm_point_labels
    )
    return kpi, band.eigenvalues
end

"""
    $(SIGNATURES)

Return the symmetry indices and labels.
"""
function get_symm_point_indices_labels(kpi::KPathInterpolant)
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
    $(SIGNATURES)

Write `prefix_band.dat, prefix_band.kpt, prefix_band.labelinfo.dat`.

This is a more user-friendly version that works with `KPathInterpolant`;
the `WannierIO.write_w90_band(prefix; kwargs...)` is the low-level version.
"""
function write_w90_band(
    prefix::AbstractString, kpi::KPathInterpolant, eigenvalues::AbstractVector
)
    kpoints = get_kpoints(kpi)
    x = get_linear_path(kpi)
    symm_point_indices, symm_point_labels = get_symm_point_indices_labels(kpi)
    return WannierIO.write_w90_band(
        prefix; x, eigenvalues, kpoints, symm_point_indices, symm_point_labels
    )
end

"""
    $(SIGNATURES)

Write kpoints into wannier90 formats: `prefix_band.kpt`, `prefix_band.labelinfo.dat`.

# Arguments
- `prefix`: the prefix of the output files `prefix_band.kpt` and `prefix_band.labelinfo.dat`
- `kpi`: a `KPathInterpolant` object

!!! tip

    This allows auto generating the high-symmetry kpoints and labels from crystal
    structure using `generate_kpath`, writing them into files. Then other codes
    can use the kpoints for band structure calculations, e.g., QE `pw.x` `bands`
    calculation, or in the `win` input file for `Wannier90`.

# Example
```julia
win = read_win("si2.win")
kp = generate_kpath(win.unit_cell_cart, win.atoms_frac, win.atom_labels)
kpi = Wannier.generate_w90_kpoint_path(kp, 100)
Wannier.write_w90_kpt_label("si2", kpi)
```
"""
function write_w90_kpt_label(prefix::AbstractString, kpi::KPathInterpolant)
    kpoints = get_kpoints(kpi)
    x = get_linear_path(kpi)
    symm_point_indices, symm_point_labels = get_symm_point_indices_labels(kpi)

    filename = prefix * "_band.kpt"
    WannierIO.write_w90_band_kpt(filename; kpoints)

    filename = prefix * "_band.labelinfo.dat"
    return WannierIO.write_w90_band_labelinfo(
        filename; x, kpoints, symm_point_indices, symm_point_labels
    )
end
