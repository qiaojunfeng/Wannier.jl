using WannierIO: SymOp, RepMatWann, RepMatBand

"""
Ensure the eigenvalue of representation matrices are integers when possible.
"""
function rescale(rep::RepMatBand)
    nbnd = size(rep.d, 1)

    d = zeros(eltype(rep.d), nbnd, nbnd)
    d .= rep.d

    for n in 1:nbnd
        s = sum(abs.(rep.d[n, :]))
        a = abs(rep.d[n, n])
        if (abs(s - a) < 1e-8) && (round(a) != 0)
            d[n, n] = rep.d[n, n] * round(a) / a
        end
    end
    return RepMatBand{nbnd}(rep.ik_ibz, rep.isym, d)
end

function rescale!(reps::AbstractVector{<:RepMatBand})
    reps .= rescale.(reps)
end

function rotate_kpoint(k::AbstractVector, symop::SymOp)
    Rk = symop.R * k
    if symop.time_reversal
        Rk = -Rk
    end
    return Rk
end

"""
    $(SIGNATURES)

Find the index mappings from kpoint in FBZ to kpoint in IBZ.

# Arguments
- `kpoints_fbz`: vector of fractional coordinates in full Brillouin zone.
- `kpoints_ibz`: vector of fractional coordinates in irreducible Brillouin zone.
- `symops`: vector of symmetry operations.

# Keyword Arguments
- `check`: whether to check the completeness and uniqueness of the mapping.

# Return
- `fbz2ibz`: a length-`nk_fbz` vector, where the `ik_fbz`-th element is a vector
    of `fbz2ibz[ik_fbz] = [ik_ibz, isym]` pairs, indicating that applying the
    `isym`-th symmetry operation to the `ik_ibz`-th kpoint in IBZ brings it to
    the `ik_fbz`-th kpoint in FBZ.
"""
function get_kpoint_mappings(
    kpoints_fbz::AbstractVector,
    kpoints_ibz::AbstractVector,
    symops::AbstractVector{SymOp};
    check::Bool=true,
)
    # For each FBZ kpoint, store a [ik_ibz, isym] pair
    fbz2ibz = [[0, 0] for _ in 1:length(kpoints_fbz)]

    # kf: FBZ kpoint
    for (ikf, kf) in enumerate(kpoints_fbz)
        # ki: IBZ kpoint
        for (iki, ki) in enumerate(kpoints_ibz)
            for (is, S) in enumerate(symops)
                Sk = rotate_kpoint(ki, S)
                if isequiv(Sk, kf)
                    j = fbz2ibz[ikf][1]
                    if j == 0
                        fbz2ibz[ikf] = [iki, is]
                    elseif j == iki
                        continue
                    else
                        error(
                            "Multiple IBZ kpoints map to the same FBZ kpoint ($ikf, $iki, $is)",
                        )
                    end
                end
            end
        end
    end

    check || return fbz2ibz

    idxs = first.(fbz2ibz)
    if 0 in idxs
        error("Some FBZ kpoints cannot be mapped to IBZ kpoints")
    end
    if sort(unique(idxs)) != collect(1:length(kpoints_ibz))
        error("Some IBZ kpoints are not mapped from any FBZ kpoints")
    end
    return fbz2ibz
end

"""
    $(SIGNATURES)

For each kpoint, store a `equiv[ik, isym]` mapping such that
`symops[isym] * kpoints[equiv[ik][isym]] = kpoints[ik]`.

# Arguments
- `kpoints`: vector of fractional coordinates.
- `symops`: vector of symmetry operations.
"""
function get_equivalence_mappings(kpoints::AbstractVector, symops::AbstractVector{SymOp})
    equiv = zeros(Int, length(kpoints), length(symops))

    # kf: final kpoint
    for (ikf, kf) in enumerate(kpoints)
        # ki: initial kpoint
        for (iki, ki) in enumerate(kpoints)
            for (is, S) in enumerate(symops)
                Sk = rotate_kpoint(ki, S)
                # Compare strictly without periodicity
                if isapprox(Sk, kf)
                    equiv[ikf, is] = iki
                end
            end
        end
    end
    if any(iszero, equiv)
        error("Some kpoints do not have equivalent mappings under all symmetry operations")
    end
    return equiv
end

"""
    $(SIGNATURES)

Unfold eigenvalues from IBZ to FBZ.

# Arguments
- `eigvals_ibz`: vector of eigenvalues at each IBZ kpoint.
- `fbz2ibz`: output of `get_kpoint_mappings`.

# Return
- `eigvals_fbz`: vector of eigenvalues at each FBZ kpoint.
"""
function unfold_eigvals(eig_ibz::AbstractVector, fbz2ibz::AbstractVector)
    nkpts_fbz = length(fbz2ibz)
    nbands = length(eig_ibz[1])
    eig = [zeros(Float64, nbands) for _ in 1:nkpts_fbz]

    for ik in 1:nkpts_fbz
        ik_ibz = fbz2ibz[ik][1]
        eig[ik] = eig_ibz[ik_ibz]
    end

    return eig
end

"""
    $(SIGNATURES)

Find the real-space translation vectors ``R_{n^\\prime}`` for each symmetry
operation ``\\hat{g}``.

The WFs are transformed according to the ``D`` matrices (around CPC Eq. 9):
```math
\\hat{g} w_{ n \\mathbf{0} }(\\mathbf{r}) =
\\sum_{n^\\prime} D_{n^\\prime n}(\\hat{g})
w_{n^\\prime \\mathbf{0}}( \\mathbf{r} - \\mathbf{R}_{n^\\prime} )
```

# Arguments
- `centers`: vector of Wannier function centers in fractional coordinates.
- `symops`: vector of symmetry operations.
- `repmat`: vector of representation matrices acting on Wannier functions.

# Return
- `Rs`: a vector of length `n_symops`, each element is a vector of length `n_wann`,
    where `Rs[is][iw]` is the translation vector for the `iw`-th Wannier function
    under the `is`-th symmetry operation.
"""
function find_wf_symmetry_translations(
    centers::AbstractVector,
    symops::AbstractVector{SymOp},
    repmat::AbstractVector{<:RepMatWann},
)
    nwann = length(centers)
    nsymm = length(symops)

    # Find where D rotates the Wannier centers
    # The translation vectors of Wannier centers, rotated WF - original WF
    Rs = [[zeros(Int, 3) for _ in 1:nwann] for _ in 1:nsymm]
    for is in 1:nsymm
        for iw in 1:nwann
            # Let's say the symmetry operation is ĝ = {S|t},
            # the Wannier centers as w,
            # we calculate rotated center S * w .
            #
            # First, find indices of the starting WFs that have connections
            # to the target iw-th WF
            jws = findall(!iszero, repmat[is].D[:, iw])
            # These WFs should have the same centers
            c = centers[jws[1]]
            @assert all(isapprox(c), [centers[jw] for jw in jws])
            # Then compute the rotated center
            # TODO check if QE source code is this convention
            # Note that the transformation on a vector r, is
            # ĝ r = r S - t
            Sw = symops[is].R' * centers[iw] - symops[is].t  # TODO very strange
            # Sw = symops[is].R * (centers[iw] + symops[is].t)
            d = Sw - c
            if all(isinteger.(d))
                Rs[is][iw] = round.(d)
            else
                error(
                    "Cannot find integer translation vector for WF $iw under symmetry $is"
                )
            end
        end
    end
    return Rs
end

"""
    $(SIGNATURES)

Unfold the gauge matrix at a *single* IBZ kpoint to the corresponding FBZ kpoint.

CPC Eq. 9.
```math
U_{m n k_f}
= \\braket{\\psi_{m k_i} | g_{0}^{-1}(k_f) | g_n}
= \\sum_j U_{m j k_i} e^{-i k_i R_j} D_{j n}(g_{0}^{-1})
```

# Arguments
- `Ui`: gauge matrix at IBZ kpoint.
- `ki`: fractional coordinate of the IBZ kpoint.
- `D`: representation matrix acting on the Wannier functions.
- `R`: vector of translation vectors for each Wannier function, fractional coordinates (integers).
- `time_reversal`: whether the symmetry operation involves time-reversal.
"""
function unfold_gauge(
    Ui::AbstractMatrix,
    ki::AbstractVector,
    D::AbstractMatrix,
    R::AbstractVector,
    time_reversal::Bool,
)
    # These all = n_wann
    size(Ui, 2) == size(D, 1) == size(D, 2) == length(R) ||
        error("Mismatch in size of Ui, D, R")

    phases = [exp(-im * 2π * dot(ki, Ri)) for Ri in R]
    # TODO multiply phases on row or column?
    # Uf = Ui * (phases .* D)
    Uf = Ui * (D .* transpose(phases))
    if time_reversal
        Uf = conj.(Uf)
    end
    return Uf
end

"""
    $(SIGNATURES)

Unfold *all* the gauge matrices from IBZ to FBZ.

CPC Eq. 9.
```math
U_{m n k_f}
= \\braket{\\psi_{m k_i} | g_{0}^{-1}(k_f) | g_n}
= \\sum_j U_{m j k_i} e^{-i k_i R_j} D_{j n}(g_{0}^{-1})
```

"""
function unfold_gauges(
    U_ibz::AbstractVector{<:AbstractMatrix},
    kpoints_ibz::AbstractVector,
    fbz2ibz::AbstractVector,
    symops::AbstractVector{SymOp},
    repmat_wann::AbstractVector{<:RepMatWann},
    Rs::AbstractVector,
)
    nk_fbz = length(fbz2ibz)
    nband, nwann = size(U_ibz[1])
    U_fbz = zeros_gauge(ComplexF64, nk_fbz, nband, nwann)

    for ik in 1:nk_fbz
        # `is` moves ik_ibz to ik_fbz
        ik_ibz, is = fbz2ibz[ik]
        ki = kpoints_ibz[ik_ibz]
        # In CPC Eq. 9, g_0(k_f) moves k_i to k_f, therefore,
        # we need the index for its inverse g_0^{-1}(k_f).
        # TODO check this inverse needed or not
        # isinv = symops[is].isym_inv
        # For now I keep it the same as original
        isinv = is
        #
        D = repmat_wann[isinv].D
        R = Rs[isinv]
        t_rev = symops[isinv].time_reversal
        U_fbz[ik] = unfold_gauge(U_ibz[ik_ibz], ki, D, R, t_rev)
    end
    return U_fbz
end

"""
    $(SIGNATURES)

Symmetrize the gauge matrices.

CPC Eq. 13.

```math
U_{m n k} = \\braket{ \\psi_{m k} | g_n }
= \\frac{1}{N_h} \\sum_h \\braket{\\psi_{m k} | h h^{-1} | g_n}
= \\frac{1}{N_h} \\sum_h \\braket{\\psi_{m k} | h | \\psi_{l k}}
                         \\braket{\\psi_{l k} | h^{-1} | g_n}
= \\frac{1}{N_h} \\sum_h d_{m l}(h, k) \\sum_j U_{l j k} e^{-i k R_j} D_{j n}(h^{-1})
```
Note that ``h \\in G_k``, the little group of kpoint.

# Return
- `U_sym`: symmetrized gauge matrices at each IBZ kpoint.
"""
function symmetrize_gauges(
    U_ibz::AbstractVector,
    kpoints_ibz::AbstractVector,
    symops::AbstractVector{SymOp},
    repmat_band::AbstractVector{<:RepMatBand},
    repmat_wann::AbstractVector{<:RepMatWann},
    Rs::AbstractVector,
)
    nband = size(repmat_band[1].d, 1)
    nwann = size(repmat_wann[1].D, 1)
    nk_ibz = length(kpoints_ibz)
    nsym = length(symops)

    U_sym = zeros_gauge(ComplexF64, nk_ibz, nband, nwann)

    idx_repmat_band = WannierIO.build_mapping_ik_isym(
        repmat_band; nkpts_ibz=nk_ibz, n_symops=nsym
    )

    for ik in 1:nk_ibz
        nh = 0
        for is in 1:nsym
            # TODO check do we need inverse here
            # isinv = symops[is].isym_inv
            # for now I keep it as original
            isinv = is
            #
            ih = idx_repmat_band[ik][is]
            isnothing(ih) && continue
            nh += 1

            Uf = unfold_gauge(
                U_ibz[ik],
                kpoints_ibz[ik],
                repmat_wann[isinv].D,
                Rs[isinv],
                symops[isinv].time_reversal,
            )
            U_sym[ik] += repmat_band[ih].d * Uf
        end
        U_sym[ik] ./= nh
    end

    return U_sym
end

"""
    $(SIGNATURES)

Get the equivalent symmetry operation of the merged operations.

This is used for CPC Eq. 19.
```math
\\hat{h} = \\hat{g}_0^{-1}(k_i + b_i) * \\hat{g}_0^{-1}(k_f) * \\hat{g}(k_f + b_f)
```

# Arguments
- `spinors`: whether spinors symmetry operations are used.
- `symops`: vector of all the symmetry operations.
- `ops`: list of operations to be merged, each element is a `isym`,
    the index of symmetry operation in `symops`.
- `invs`: list of whether to take inverse of each operation, each element is a bool,
    `true` means take inverse.

# Return
- `isym_h`: index of the equivalent symmetry operation in `symops`.
- `T`: translation vector, the difference between the l.h.s and r.h.s. of the
    previous equation.
"""
function merge_symops(
    spinors::Bool,
    symops::AbstractVector{SymOp},
    ops::AbstractVector{<:Integer},
    invs::AbstractVector{Bool},
)
    # Initialize everything to identity
    s0 = Matrix{Float64}(I, 3, 3)
    t0 = zeros(3)
    u0 = Matrix{ComplexF64}(I, 2, 2)
    t_rev = false

    # -i * σ_y, also need conjugate
    Trev = [0 1; -1 0]
    # (-i σ_y)^-1
    Trev_inv = [0 -1; 1 0]

    for (isym, inv_op) in zip(ops, invs)
        if inv_op
            # inv(r*S-t) = (r+t)*S^-1 = r*S^-1 + t*S^-1
            s1 = symops[symops[isym].isym_inv].R
            t1 = - transpose(s1) * symops[isym].t
            u1 = symops[isym].u'
            if symops[isym].time_reversal
                u1 = u1 * Trev_inv
            end
        else
            # r*S-t
            s1 = symops[isym].R
            t1 = symops[isym].t
            u1 = symops[isym].u
            if symops[isym].time_reversal
                u1 = Trev * conj.(u1)
            end
        end

        # Now we merge operation 0 and 1
        # r' = r*S0 - t0
        # r''= r'*S1 - t1 = (r*S0 - t0)*S1 - t1 = r*S0*S1 - t0*S1 - t1
        s0 *= s1
        t0 = transpose(s1) * t0 + t1
        if t_rev
            u1 = conj.(u1)
        end
        u0 *= u1
        t_rev = xor(t_rev, symops[isym].time_reversal)
    end

    # Find operation in symops
    for (isym, op) in enumerate(symops)
        T = t0 - op.t
        if isapprox(s0, op.R; atol=1e-6) &&
            isapprox(T, round.(T); atol=1e-6) &&
            (t_rev == op.time_reversal)
            T = Int.(T)
            if spinors
                if t_rev
                    us = Trev * conj.(symops[isym].u)
                else
                    us = symops[isym].u
                end
                if isapprox(u0, us; atol=1e-5)
                    return isym, 1, T
                elseif isapprox(u0, -us; atol=1e-5)
                    return isym, -1, T
                else
                    error("u does not match")
                end
            else
                return isym, 1, T
            end
        end
    end
    return error("No equivalent symmetry operation found")
end

"""
    $(SIGNATURES)

Unfold overlap matrices from IBZ to FBZ.

CPC Eq. 6 and 7.

```math
M_{m n}^{k_f, b_f} = \\sum_l M_{m l}^{k_i, b_i} d_{l n}(\\hat{h}, k_i)
```

# Arguments
- `M_ibz`: overlap matrices at each IBZ kpoint.
- `kpoints_ibz`: fractional coordinates of IBZ kpoints.
- `fbz2ibz`: output of `get_kpoint_mappings`.
- `kstencil`: k-space stencil for the FBZ.
- `spinors`: whether spinors symmetry operations are used.
- `symops`: vector of symmetry operations.
- `repmat_band`: representation matrices acting on the Bloch states.

# Return
- `M_fbz`: overlap matrices at each FBZ kpoint. The b vectors are ordered
    according to `kstencil`.
"""
function unfold_overlaps(
    M_ibz::AbstractVector,
    kpoints_ibz::AbstractVector,
    fbz2ibz::AbstractVector,
    kstencil::KspaceStencil,
    spinors::Bool,
    symops::AbstractVector{SymOp},
    repmat_band::AbstractVector{<:RepMatBand},
)
    nk_fbz = n_kpoints(kstencil)
    length(fbz2ibz) == nk_fbz || error("Mismatch in number of FBZ kpoints")
    nband = size(M_ibz[1][1], 1)
    nbvec = n_bvectors(kstencil)
    kpoints_fbz = kstencil.kpoints
    Mf = zeros_overlap(ComplexF64, nk_fbz, nbvec, nband)

    ikisym2ih = WannierIO.build_mapping_ik_isym(
        repmat_band; nkpts_ibz=length(kpoints_ibz), n_symops=length(symops)
    )

    for (ikf, kf) in enumerate(kpoints_fbz)
        # isym: from ki to kf
        iki, isym_kf = fbz2ibz[ikf]
        ki = kpoints_ibz[iki]
        bvecs = get_bvectors(kstencil, ikf; fractional=true)
        # Get mapping between equivalent b vectors
        b2b = get_equivalence_mappings(bvecs, symops)

        for (ibf, bf) in enumerate(bvecs)
            # In CPC Eq. 6, g_0(k_f) is from ki to kf, we need its inverse
            # g_0^{-1}(k_f) to go from kf to ki, to obtain b_i
            # b_i = R^{-1} b_f where g_0 = {R|t}
            # bi = rotate_kpoint(bf, symops[symops[isym_kf].isym_inv])
            # ibi = index_bvector(kpoints_ibz, kpb_k_ibz, kpb_G_ibz, iki, bi)
            # isnothing(ibi) && error("No bvector found for ik = $(iki), bi = $(bi)")
            ibi = b2b[ibf, isym_kf]
            bi = bvecs[ibi]

            # We need symmetries at ki + bi:
            # 1. find the index of ki + bi in the FBZ
            ikbi_fbz = findfirst(isequiv(ki + bi), kpoints_fbz)
            # 2. find the index of ki+bi in the IBZ, and symmetry to get ki+bi in IBZ
            ikbi_ibz, isym_kbi = fbz2ibz[ikbi_fbz]
            # 3. find the index of kf + bf in the FBZ
            ikbf_fbz = findfirst(isequiv(kf + bf), kpoints_fbz)
            # 4. find the index of kf+bf in the IBZ, and symmetry to get kf+bf in IBZ
            isym_kbf = fbz2ibz[ikbf_fbz][2]

            # get equivalent operation of h = g₀⁻¹(ki+bi) * g₀⁻¹(kf) * g₀(kf+bf)
            isym_h, factor, T = merge_symops(
                spinors, symops, [isym_kbi, isym_kf, isym_kbf], [true, true, false]
            )

            ih = ikisym2ih[ikbi_ibz][isym_h]
            d = repmat_band[ih].d
            if symops[isym_kbi].time_reversal
                Mf[ikf][ibf] = M_ibz[iki][ibi] * conj.(d) * factor
            else
                Mf[ikf][ibf] = M_ibz[iki][ibi] * d * factor
            end
            if symops[isym_kf].time_reversal
                Mf[ikf][ibf] = conj.(Mf[ikf][ibf])
            end

            kbi_ibz = kpoints_ibz[ikbi_ibz]
            θ1 = dot(bi, symops[isym_kf].t)
            θ2 = dot(kbi_ibz, T)
            if symops[isym_kbi].time_reversal
                θ2 *= -1
            end
            if symops[isym_kf].time_reversal
                θ1 *= -1
                θ2 *= -1
            end
            if symops[isym_h].time_reversal
                θ2 *= -1
            end
            phase = exp(-im * 2π * (θ1 + θ2))
            Mf[ikf][ibf] *= phase
        end
    end
    return Mf
end
