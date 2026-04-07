export unfold

"""
    $(SIGNATURES)

Unfold the files computed in the irreducible Brillouin zone (IBZ) to the
full Brillouin zone (FBZ).

Ref: T. Koretsune, CPC 2023, 108645.
https://doi.org/10.1016/j.cpc.2022.108645

# Arguments
- `prefix`: prefix of the input files, e.g. `wannier`.
- `out_prefix`: prefix of the output files, default to be the same as `prefix`.

# Keyword Arguments
- `reorder_bvec`: Whether to reorder the b vectors such that the output
    `mmn` file has the same ordering as the input `nnkp` file.
    Default is `false`, i.e., the b vectors are always the same across all the
    kpoints (same as the Γ point).

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
function unfold(
        prefix::AbstractString, out_prefix::AbstractString = prefix; reorder_bvec::Bool = false
    )
    nnkp = read_nnkp("$prefix.nnkp")
    kstencil = Wannier.KspaceStencil(
        nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
    )

    isym = read_isym("$prefix.isym")
    # TODO really needed?
    Wannier.rescale!(isym.repmat_band)

    kpoints_ibz = isym.kpoints_ibz
    symops = isym.symops
    repmat_band = isym.repmat_band
    repmat_wann = isym.repmat_wann

    f2i = get_kpoint_mappings(kstencil.kpoints, kpoints_ibz, symops)

    # eig
    Ei = read_eig("$prefix.ieig")
    Ef = Wannier.unfold_eigvals(Ei, f2i)
    write_eig("$out_prefix.eig", Ef)

    # amn
    Ai = WannierIO.read_amn("$prefix.iamn").A
    # The factor exp(-i kᵢ R_{n'}) appearing in CPC Eq. 9
    centers = [p.center for p in nnkp["projections"]]
    Rs = Wannier.find_wf_symmetry_translations(centers, symops, repmat_wann)
    Asymm = Wannier.symmetrize_gauges(Ai, kpoints_ibz, symops, repmat_band, repmat_wann, Rs)
    Af = Wannier.unfold_gauges(Asymm, kpoints_ibz, f2i, symops, repmat_wann, Rs)
    write_amn("$out_prefix.amn", Af)

    # mmn
    mmn_i = read_mmn("$prefix.immn")
    Mi = mmn_i.M
    kpb_k_i = mmn_i.kpb_k
    kpb_G_i = mmn_i.kpb_G
    # The b vectors of the IBZ mmn file is the same as that of the FBZ nnkp file
    # at the Γ point, and the b vectors of remaining kpoints are always the same
    # at each kpoint in the IBZ mmn file.
    kstencil_ibz = Wannier.KspaceStencil(
        kstencil.recip_lattice, kstencil.kgrid_size, kpoints_ibz,
        get_bvectors(kstencil; fractional = true),
        kstencil.bweights, kpb_k_i, kpb_G_i
    )
    Mf, kpb_k_f, kpb_G_f = Wannier.unfold_overlaps(
        Mi, kstencil_ibz, kstencil, f2i, isym.spinors, symops, repmat_band
    )

    if reorder_bvec
        Mf = Wannier.reorder(Mf, kpb_k_f, kpb_G_f, kstencil)
    else
        kstencil.kpb_k .= kpb_k_f
        kstencil.kpb_G .= kpb_G_f
    end
    write_mmn("$out_prefix.mmn", Mf, kstencil)

    return nothing
end
