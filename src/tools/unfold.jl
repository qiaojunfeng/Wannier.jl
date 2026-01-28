export unfold

"""
    $(SIGNATURES)

Unfold the files computed in the irreducible Brillouin zone (IBZ) to the
full Brillouin zone (FBZ).

Ref: T. Koretsune, CPC 2023, 108645.
https://doi.org/10.1016/j.cpc.2022.108645

# Arguments
- `prefix`: prefix of the input files, e.g. `wannier`.

# Input files
- `prefix.nnkp`: nearest neighbor kpoint file.
- `prefix.isym`: symmetry operation file.
- `prefix.ieig`: eigenvalue file in IBZ.
- `prefix.iamn`: projection matrix in IBZ.
- `prefix.immn`: overlap matrix element file in IBZ.

# Output
- `prefix.eig`: eigenvalue file in FBZ.
- `prefix.amn`: projection matrix in FBZ.
- `prefix.mmn`: overlap matrix in FBZ.
"""
function unfold(prefix::AbstractString)
    nnkp = read_nnkp("$prefix.nnkp")
    # kstencil = read_nnkp_compute_bweights("$prefix.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G
    )

    isym = read_isym("$prefix.isym")
    # TODO really needed?
    isym.repmat_band = rescale.(isym.repmat_band)

    kpoints_ibz = isym.kpoints_ibz
    symops = isym.symops
    repmat_band = isym.repmat_band
    repmat_wann = isym.repmat_wann

    f2i = get_kpoint_mappings(kstencil.kpoints, kpoints_ibz, symops)

    # eig
    Ei = read_eig("$prefix.ieig")
    Ef = unfold_eigvals(Ei, f2i)
    write_eig("$prefix.eig", Ef)

    # amn
    Ai = read_amn("$prefix.iamn")
    # The factor exp(-i kᵢ R_{n'}) appearing in CPC Eq. 9
    centers = [p.center for p in nnkp.projections]
    Rs = find_wf_symmetry_translations(centers, symops, repmat_wann)
    Asymm = symmetrize_gauges(Ai, kpoints_ibz, symops, repmat_band, repmat_wann, Rs)
    Af = unfold_gauges(Asymm, kpoints_ibz, f2i, symops, repmat_wann, Rs)
    write_amn("$prefix.amn", Af)

    # mmn
    Mi, kpb_k_i, kpb_G_i = read_mmn("$prefix.immn")
    Mf = unfold_overlaps(
        Mi, kpb_k_i, kpb_G_i, kpoints_ibz, f2i, kstencil, symops, repmat_band
    )
    return write_mmn("$prefix.mmn", Mf, kstencil)
end
