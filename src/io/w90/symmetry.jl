export read_w90_ibz

"""
    read_w90_ibz(prefix; centers=nothing, frozen_band_indices=nothing,
                 frozen_bands=nothing, orthonormalize_projections=true,
                 clean_representations=true)

Read the irreducible-mesh Wannier90 interface files and construct a
[`SymmetricModel`](@ref) ready for symmetry-constrained localization.

The files `prefix.nnkp`, `prefix.isym`, `prefix.ieig`, `prefix.immn`,
`prefix.iamn`, and `prefix.win` must exist. The reader normalizes the
little-group matrices, constructs the symmetry constraint, reconstructs the
full-mesh overlaps, eigenvalues, and initial gauge, and uses the global
b-vector ordering required by the symmetry kernels.

By default, WF centers come from the projections in `.nnkp`, and frozen bands
come from `dis_froz_min`/`dis_froz_max` in `.win`. Pass exact fractional
`centers` when the finite precision of `.nnkp` does not preserve a special
Wyckoff coordinate. A fixed band selection can be supplied through
`frozen_band_indices`; alternatively, `frozen_bands` accepts a complete
`n_bands × n_kpoints_fbz` Boolean mask.
"""
function read_w90_ibz(
        prefix::AbstractString;
        centers = nothing,
        frozen_band_indices = nothing,
        frozen_bands = nothing,
        orthonormalize_projections::Bool = true,
        clean_representations::Bool = true,
    )
    required = ("nnkp", "isym", "ieig", "immn", "iamn", "win")
    for extension in required
        path = "$prefix.$extension"
        isfile(path) || error("$path file does not exist")
    end

    nnkp = read_nnkp("$prefix.nnkp")
    stencil_ibz = KSpaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )
    isym = read_isym("$prefix.isym")
    eigenvalues_ibz = read_eig("$prefix.ieig")
    if clean_representations
        clean_littlegroup_reps!(isym.littlegroup_reps, eigenvalues_ibz)
    else
        normalize_diagonal_littlegroup_reps!(isym.littlegroup_reps)
    end

    target_centers = if isnothing(centers)
        [projection.center for projection in nnkp["projections"]]
    else
        collect(centers)
    end
    constraint = SymmetryConstraint(stencil_ibz, isym, target_centers)

    overlaps_ibz = read_mmn("$prefix.immn").M
    overlaps_fbz = reconstruct_overlaps(overlaps_ibz, constraint)
    eigenvalues_fbz = unfold_eigvals(
        eigenvalues_ibz, [collect(pair) for pair in constraint.fbz2ibz]
    )

    projections_ibz = read_amn("$prefix.iamn").A
    projections_ibz = project_covariant(projections_ibz, constraint)
    projections_fbz = reconstruct_gauges(projections_ibz, constraint)
    if orthonormalize_projections
        projections_fbz = lowdin_orthonormalize(projections_fbz)
    end

    win = read_win("$prefix.win")
    nbands, nkpoints_fbz = size(eigenvalues_fbz)
    frozen = if !isnothing(frozen_bands)
        isnothing(frozen_band_indices) || throw(
            ArgumentError("pass either frozen_bands or frozen_band_indices, not both")
        )
        size(frozen_bands) == size(eigenvalues_fbz) || throw(
            DimensionMismatch(
                "frozen_bands must have size $(size(eigenvalues_fbz)), got $(size(frozen_bands))"
            )
        )
        BitMatrix(frozen_bands)
    elseif !isnothing(frozen_band_indices)
        frozen_mask = falses(nbands, nkpoints_fbz)
        frozen_mask[frozen_band_indices, :] .= true
        frozen_mask
    else
        frozen_max = get(win, "dis_froz_max", nothing)
        if isnothing(frozen_max)
            falses(nbands, nkpoints_fbz)
        else
            get_frozen_bands(
                eigenvalues_fbz, frozen_max, get(win, "dis_froz_min", -Inf)
            )
        end
    end

    atom_positions = [atom.second for atom in win["atoms_frac"]]
    atom_labels = [string(atom.first) for atom in win["atoms_frac"]]
    model = Model(
        win["unit_cell_cart"],
        atom_positions,
        atom_labels,
        globalize_bvector_ordering(stencil_ibz),
        overlaps_fbz,
        projections_fbz,
        eigenvalues_fbz,
        frozen,
    )
    return SymmetricModel(model, constraint, overlaps_ibz)
end
