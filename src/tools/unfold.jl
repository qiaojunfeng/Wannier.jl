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
    kstencil = KspaceStencil(nnkp.recip_lattice, nnkp.kpoints, nnkp.kpb_k, nnkp.kpb_G)

    isym = read_isym("$prefix.isym")
    # TODO really needed?
    rescale!.(isym.repmat_band)

    f2i = get_kpoint_mappings(kstencil.kpoints, isym.kpoints_ibz, isym.symops)

    # eig
    Ei = read_eig("$prefix.ieig")
    Ef = unfold_eigvals(Ei, f2i)

    # amn
    Ai = read_amn("$prefix.iamn")
    # The factor exp(-i kᵢ R_{n'}) appearing in CPC Eq. 9
    centers = [p.center for p in nnkp.projections]
    Rs = find_wf_symmetry_translations(centers, isym.symops, isym.repmat_wann)
    Asymm = symmetrize_gauges(Ai, isym.kpoints_ibz, isym.symops, isym.repmat_band, isym.repmat_wann, Rs)
    Af = unfold_gauges(Asymm, isym.kpoints_ibz, f2i, isym.symops, isym.repmat_wann, Rs)

    # mmn
    Mi, kpb_k_i, kpb_G_i = read_mmn("$prefix.immn")
    Mf = unfold_overlaps(
        Mi, kpb_k_i, kpb_G_i, isym.kpoints_ibz, f2i, kstencil, isym.symops, isym.repmat_band
    )

    #check Eig with symwannier
    nk_fbz = n_kpoints(kstencil)
    eigpy = read_eig("$prefix.eig")
    for ik in 1:nk_fbz
        i = findall(>(1e-10), abs.(eigpy[ik]-Ef[ik]))
        println("ik = ", ik, "    ", i)
    end

    #check Amn with symwannier
    Apy = read_amn("$prefix.amn")
    for ik in 1:nk_fbz
        i = findall(>(1e-10), abs.(Af[ik]-Apy[ik]))
        println("ik = ", ik, "    ", i)
    end

    #check Mmn with symwannier
    Mpy = read_mmn("$prefix.mmn")
    for ik in 1:nk_fbz
        for ib in 1:8
            i = findall(>(1e-10), abs.(Mf[ik][ib]-Mpy[ik][ib]))
            println("ik = ", ik, "    ", i)
        end
    end
end
