using LinearAlgebra
using WannierIO: SymOp, OrbitalRep, LittleGroupRep

# All symmetry operations are in the standard (ITA) Seitz convention as
# returned by `WannierIO.read_isym` (see `WannierIO.standardize`):
# - real space (fractional): ĝ r = W r + v
# - kpoints (fractional): k′ = Wk k, with Wk = transpose(inv(W));
#   for antiunitary (time-reversal) operations k′ = -Wk k
# - `littlegroup_reps` store d(ĥ, k) = ⟨ψ_m|ĥ ψ_n⟩ (column = original state)
# - `orbital_reps[isym]` store D(ĝ_isym): ĝ w_n = Σ_{n′} D_{n′n} w_{n′}(r - R_{n′})

"""
Ensure the eigenvalue of representation matrices are integers when possible.
"""
function rescale(rep::LittleGroupRep)
    nbnd = size(rep.d, 1)

    d = zeros(eltype(rep.d), nbnd, nbnd)
    d .= rep.d

    for n in 1:nbnd
        s = sum(abs.(rep.d[n, :]))
        a = abs(rep.d[n, n])
        if (abs(s - a) < 1.0e-8) && (round(a) != 0)
            d[n, n] = rep.d[n, n] * round(a) / a
        end
    end
    return LittleGroupRep{nbnd}(rep.ik_ibz, rep.isym, d)
end

function rescale_littlegroup_reps!(reps::AbstractVector{<:LittleGroupRep})
    return reps .= rescale.(reps)
end

export clean_littlegroup_reps!, map_fbz_to_ibz

"""
    $(SIGNATURES)

Remove the numerical noise from little-group representation matrices.

For energy multiplets fully inside a symmetry-closed window the matrices
`d(ĥ, k)` are *exactly* block-diagonal over degenerate multiplets (`ĥ`
commutes with the Hamiltonian) and *exactly* unitary within each block; the
`.isym` data carries them only to the accuracy of the generating DFT run
(~1e-5 entry noise is common), and that noise floors every downstream
symmetry tolerance. This function enforces both exact properties: the bands
at each rep's IBZ kpoint are clustered by energy gaps larger than
`atol_degeneracy` (the same rule as the energy-multiplet masking of
[`SymmetryConstraint`](@ref)), all cross-cluster entries are zeroed, and
each within-cluster block is replaced by its closest unitary (polar factor,
`lowdin_orthonormalize`).

A block is unitarized only when it is already unitary to within
`atol_unitary` (measured as `opnorm(B'B - I)`). Blocks with a larger deficit
belong to multiplets the window genuinely truncates (or that the data
genuinely breaks): their `d` entries are contractions, not noisy unitaries,
and they are left untouched so the symmetry-broken-band masking of
[`SymmetryConstraint`](@ref) still detects them.

`eig_ibz` are the IBZ eigenvalues (`n_bands × n_kpoints_ibz`, the `.ieig`
data, ascending per kpoint). Apply after `rescale_littlegroup_reps!`. Cleaning is
strictly opt-in: it moves quantities unfolded through the reps (e.g.
`reconstruct_overlaps`) by the size of the removed noise, so data cleaned
here no longer reproduces reference files generated with the raw reps to
better than that noise.
"""
function clean_littlegroup_reps!(
        reps::AbstractVector{<:LittleGroupRep},
        eig_ibz::AbstractMatrix{<:Real};
        atol_degeneracy::Real = 1.0e-4,
        atol_unitary::Real = 0.01,
    )
    isempty(reps) && return reps
    nbnd = size(reps[1].d, 1)
    size(eig_ibz, 1) == nbnd ||
        error("eig_ibz must have n_bands = $nbnd rows, got $(size(eig_ibz, 1))")
    return reps .= map(reps) do rep
        E = view(eig_ibz, :, rep.ik_ibz)
        d0 = Matrix(rep.d)
        d = zeros(eltype(d0), nbnd, nbnd)
        lo = 1
        for n in 1:nbnd
            if n == nbnd || E[n + 1] - E[n] > atol_degeneracy
                blk = lo:n
                B = d0[blk, blk]
                if opnorm(B' * B - I) <= atol_unitary
                    B = lowdin_orthonormalize(B)
                end
                d[blk, blk] = B
                lo = n + 1
            end
        end
        LittleGroupRep{nbnd}(rep.ik_ibz, rep.isym, d)
    end
end

function rotate_kpoint(k::AbstractVector, symop::SymOp)
    Rk = symop.Wk * k
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
function map_fbz_to_ibz(
        kpoints_fbz::AbstractVector,
        kpoints_ibz::AbstractVector,
        symops::AbstractVector{SymOp};
        check::Bool = true,
    )
    # For each FBZ kpoint, store a [ik_ibz, isym] pair
    fbz2ibz = [[0, 0] for _ in 1:length(kpoints_fbz)]

    # kf: FBZ kpoint
    for (ikf, kf) in enumerate(kpoints_fbz)
        # ki: IBZ kpoint
        for (iki, ki) in enumerate(kpoints_ibz)
            for (is, S) in enumerate(symops)
                Sk = rotate_kpoint(ki, S)
                if isequivalent(Sk, kf)
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
- `fbz2ibz`: output of `map_fbz_to_ibz`.

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

function unfold_eigvals(eig_ibz::AbstractMatrix, fbz2ibz::AbstractVector)
    nbands = size(eig_ibz, 1)
    nkpts_fbz = length(fbz2ibz)
    eig = zeros(eltype(eig_ibz), nbands, nkpts_fbz)
    for ik in 1:nkpts_fbz
        ik_ibz = fbz2ibz[ik][1]
        eig[:, ik] .= view(eig_ibz, :, ik_ibz)
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
so that ``\\hat{g} \\tau_n = \\tau_{n^\\prime} + R_{n^\\prime}(\\hat{g})``,
where ``\\tau_n`` are the WF centers.

# Arguments
- `centers`: vector of Wannier function centers in fractional coordinates.
- `symops`: vector of symmetry operations.
- `orbital_reps`: representation matrices `D(ĝ)` acting on the Wannier functions.

# Return
- `Rs`: a vector of length `n_symops`, each element is a vector of length `n_wann`,
    where `Rs[is][iw]` is the translation vector
    ``R_{n^\\prime}(\\hat{g}_{is}) = \\hat{g}_{is} \\tau_{iw} - \\tau_{n^\\prime}``
    for the row block ``n^\\prime`` connected to the `iw`-th *column* of
    ``D(\\hat{g}_{is})``.
"""
function find_wf_symmetry_translations(
        centers::AbstractVector,
        symops::AbstractVector{SymOp},
        orbital_reps::AbstractVector{<:OrbitalRep},
    )
    nwann = length(centers)
    nsymm = length(symops)

    Rs = [[zeros(Int, 3) for _ in 1:nwann] for _ in 1:nsymm]
    for is in 1:nsymm
        for iw in 1:nwann
            # Find the rows n′ connected to column iw: D_{n′,iw}(ĝ) ≠ 0 iff
            # ĝ τ_iw = τ_{n′} + R with R an integer lattice vector
            jws = findall(!iszero, orbital_reps[is].D[:, iw])
            # These WFs should have the same centers
            c = centers[jws[1]]
            @assert all(isapprox(c), [centers[jw] for jw in jws])
            # Rotated center in the standard Seitz convention: ĝ τ = W τ + v
            Sw = symops[is].W * centers[iw] + symops[is].v
            d = Sw - c
            # They should be integer translations, but `all(isinteger.(d))` is too strict
            if isapprox(d, round.(d); atol = 1.0e-8)
                Rs[is][iw] = round.(d)
            else
                error(
                    "Cannot find integer translation vector for WF $iw under symmetry $is, d = $d",
                )
            end
        end
    end
    return Rs
end

"""
    $(SIGNATURES)

Lattice translation ``L`` such that
``\\hat{g}_{isym\\_inv} \\circ \\hat{g}_{isym} = t_L``, i.e. the *stored*
inverse element differs from the exact inverse by a lattice translation:
``\\hat{g}_{isym}^{-1} = t_{-L} \\circ \\hat{g}_{isym\\_inv}``.

For symmorphic groups ``L = 0``; for nonsymmorphic groups it is generally
nonzero and must be accounted for whenever the *exact* inverse operator is
required (e.g. the translation vectors ``R_{n^\\prime}(g_0^{-1})`` in the
gauge unfolding).
"""
function inverse_translation_mismatch(symops::AbstractVector{SymOp}, isym::Integer)
    op = symops[isym]
    opinv = symops[op.isym_inv]
    L = opinv.v + opinv.W * op.v
    Li = round.(Int, L)
    isapprox(L, Li; atol = 1.0e-6) ||
        error("stored inverse is not the inverse modulo a lattice translation for isym = $isym")
    return Vec3{Int}(Li)
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
- `D`: representation matrix ``D(g_0^{-1})`` acting on the Wannier functions.
- `R`: vector of translation vectors ``R_{n^\\prime}(g_0^{-1})``, fractional
    coordinates (integers), indexed by the *column* of `D`.
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
    # `R[j]` is indexed by the *column* j of D; since D's nonzero pattern is a
    # permutation of site blocks (row block n' fixed by column j), scaling
    # column j by exp(-i k R[j]) is identical to the row scaling
    # exp(-i k R_n') in the unfolding formula.
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
        U_ibz::AbstractArray{<:Complex, 3},
        kpoints_ibz::AbstractVector,
        fbz2ibz::AbstractVector,
        symops::AbstractVector{SymOp},
        orbital_reps::AbstractVector{<:OrbitalRep},
        Rs::AbstractVector,
    )
    nk_fbz = length(fbz2ibz)
    nband, nwann = size(U_ibz, 1), size(U_ibz, 2)
    U_fbz = zeros_gauge(ComplexF64, nk_fbz, nband, nwann)

    for ik in 1:nk_fbz
        # `is` moves ik_ibz to ik_fbz, i.e. g₀(k_f) = symops[is]; the formula
        # needs D(g₀⁻¹(k_f)) and R(g₀⁻¹(k_f)), so index by the inverse.
        # The stored inverse element equals the exact inverse only up to a
        # lattice translation, g₀⁻¹ = t₋L ∘ ĝ_isinv; D is unaffected but the
        # translation vectors must be corrected: R(g₀⁻¹) = R(ĝ_isinv) - L.
        ik_ibz, is = fbz2ibz[ik]
        ki = kpoints_ibz[ik_ibz]
        isinv = symops[is].isym_inv
        D = orbital_reps[isinv].D
        L = inverse_translation_mismatch(symops, is)
        R = [Ri - L for Ri in Rs[isinv]]
        t_rev = symops[isinv].time_reversal
        view(U_fbz, :, :, ik) .= unfold_gauge(view(U_ibz, :, :, ik_ibz), ki, D, R, t_rev)
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
        U_ibz::AbstractArray{<:Complex, 3},
        kpoints_ibz::AbstractVector,
        symops::AbstractVector{SymOp},
        littlegroup_reps::AbstractVector{<:LittleGroupRep},
        orbital_reps::AbstractVector{<:OrbitalRep},
        Rs::AbstractVector,
    )
    nband = size(littlegroup_reps[1].d, 1)
    nwann = size(orbital_reps[1].D, 1)
    nk_ibz = length(kpoints_ibz)
    nsym = length(symops)

    U_sym = zeros_gauge(ComplexF64, nk_ibz, nband, nwann)

    idx_littlegroup = WannierIO.build_mapping_ik_isym(
        littlegroup_reps; nkpts_ibz = nk_ibz, n_symops = nsym
    )

    for ik in 1:nk_ibz
        nh = 0
        for is in 1:nsym
            ih = idx_littlegroup[ik][is]
            isnothing(ih) && continue
            nh += 1

            # The formula needs D(h⁻¹) and R(h⁻¹); index by the inverse, and
            # correct the translations for the stored-inverse lattice mismatch
            # (see `unfold_gauges`).
            isinv = symops[is].isym_inv
            L = inverse_translation_mismatch(symops, is)
            Uf = unfold_gauge(
                view(U_ibz, :, :, ik),
                kpoints_ibz[ik],
                orbital_reps[isinv].D,
                [Ri - L for Ri in Rs[isinv]],
                symops[isinv].time_reversal,
            )
            view(U_sym, :, :, ik) .+= littlegroup_reps[ih].d * Uf
        end
        view(U_sym, :, :, ik) ./= nh
    end

    return U_sym
end

"""
    $(SIGNATURES)

Get the equivalent symmetry operation of the merged operations.

This is used for CPC Eq. 19.
```math
\\hat{h} = \\hat{g}_0^{-1}(k_i + b_i) * \\hat{g}_0^{-1}(k_f) * \\hat{g}_0(k_f + b_f)
```
The list is composed as an operator product with the *rightmost* element
applied first, i.e. `ops = [i₁, i₂, i₃]` with `invs = [true, true, false]`
yields ``\\hat{h} = \\hat{g}_{i_1}^{-1} \\hat{g}_{i_2}^{-1} \\hat{g}_{i_3}``.

# Arguments
- `spinors`: whether spinors symmetry operations are used.
- `symops`: vector of all the symmetry operations.
- `ops`: list of operations to be merged, each element is a `isym`,
    the index of symmetry operation in `symops`.
- `invs`: list of whether to take inverse of each operation, each element is a bool,
    `true` means take inverse.

# Return
- `isym_h`: index of the equivalent symmetry operation in `symops`.
- `factor`: ±1 sign relating the composed SU(2) matrix to the stored one
    (double-group sign; always `1` for `spinors = false`).
- `T`: integer lattice translation such that the composed operation equals
    ``\\hat{g}_{isym_h} \\circ t_T``, i.e. ``r \\mapsto \\hat{g}_{isym_h}(r + T)``.
"""
function compose_symops(
        spinors::Bool,
        symops::AbstractVector{SymOp},
        ops::AbstractVector{<:Integer},
        invs::AbstractVector{Bool},
    )
    # Accumulate the operator product h = f₁ ∘ f₂ ∘ ... (rightmost applied
    # first): appending the next factor f gives acc ∘ f, i.e.
    # {W₀|v₀} ∘ {W₁|v₁} = {W₀W₁ | v₀ + W₀v₁}.
    W0 = Mat3{Int}(I)
    v0 = zeros(3)
    u0 = Matrix{ComplexF64}(I, 2, 2)
    t_rev = false

    # -i * σ_y, also need conjugate
    Trev = [0 1; -1 0]
    # (-i σ_y)^-1
    Trev_inv = [0 -1; 1 0]

    for (isym, inv_op) in zip(ops, invs)
        op = symops[isym]
        if inv_op
            # {W|v}⁻¹ = {W⁻¹ | -W⁻¹v}; W⁻¹ is the (exact, integer) rotation
            # of the stored inverse element.
            W1 = symops[op.isym_inv].W
            v1 = -(W1 * op.v)
            u1 = op.u'
            if op.time_reversal
                u1 = u1 * Trev_inv
            end
        else
            W1 = op.W
            v1 = op.v
            u1 = op.u
            if op.time_reversal
                u1 = Trev * conj.(u1)
            end
        end

        # Spatial composition acc ∘ f (θ̂ commutes with spatial operations)
        v0 = v0 + W0 * v1
        W0 = W0 * W1
        # SU(2) composition with the magnetic-group rule
        # u(acc ∘ f) = u_acc ⋅ K_acc(u_f)
        if t_rev
            u1 = conj.(u1)
        end
        u0 *= u1
        t_rev = xor(t_rev, op.time_reversal)
    end

    # Find the stored operation: the composite equals ĝ_isym ∘ t_T with
    # T = W_isym⁻¹ (v₀ - v_isym) an integer lattice vector.
    for (isym, op) in enumerate(symops)
        op.time_reversal == t_rev || continue
        W0 == op.W || continue
        T = symops[op.isym_inv].W * (v0 - op.v)
        isapprox(T, round.(T); atol = 1.0e-6) || continue
        T = Int.(round.(T))
        if spinors
            if t_rev
                us = Trev * conj.(op.u)
            else
                us = op.u
            end
            if isapprox(u0, us; atol = 1.0e-5)
                return isym, 1, T
            elseif isapprox(u0, -us; atol = 1.0e-5)
                return isym, -1, T
            else
                error("u does not match")
            end
        else
            return isym, 1, T
        end
    end
    return error("No equivalent symmetry operation found")
end

"""
    $(SIGNATURES)

Reconstruct overlap matrices from the IBZ on the full Brillouin-zone mesh.

CPC Eq. 6 and 7.

```math
M_{m n}^{k_f, b_f} = \\sum_l M_{m l}^{k_i, b_i} d_{l n}(\\hat{h}, k_i)
```

# Arguments
- `M_ibz`: overlap matrices at each IBZ kpoint.
- `kpb_k_ibz`: index of k+b at each IBZ kpoint.
- `kpb_G_ibz`: G-vector part of k+b at each IBZ kpoint.
- `kpoints_ibz`: fractional coordinates of IBZ kpoints.
- `bvectors`: b vectors in fractional coordinates.
    The IBZ mmn has the same b vectors ordering for all the kpoints.
- `kpoints_fbz`: fractional coordinates of FBZ kpoints.
- `fbz2ibz`: output of `map_fbz_to_ibz`.
- `spinors`: whether spinors symmetry operations are used.
- `symops`: vector of symmetry operations.
- `littlegroup_reps`: representation matrices acting on the Bloch states.

# Return
- `M_fbz`: overlap matrices at each FBZ kpoint. The b vectors are ordered
    according to `kstencil`.
"""
function reconstruct_overlaps(
        M_ibz::AbstractArray{<:Complex, 4},
        kpb_k_ibz::AbstractMatrix{Int},
        kpb_G_ibz::AbstractMatrix,
        kpoints_ibz::AbstractVector,
        bvectors::AbstractVector,
        kpoints_fbz::AbstractVector,
        fbz2ibz::AbstractVector,
        spinors::Bool,
        symops::AbstractVector{SymOp},
        littlegroup_reps::AbstractVector{<:LittleGroupRep},
    )
    nk_fbz = length(kpoints_fbz)
    length(fbz2ibz) == nk_fbz || error("Mismatch in number of FBZ kpoints")
    nband = size(M_ibz, 1)
    nbvec = length(bvectors)
    Mf = zeros_overlap(ComplexF64, nk_fbz, nbvec, nband)
    kpb_k_fbz = zeros(Int, nbvec, nk_fbz)
    kpb_G_fbz = Matrix{Vec3{Int}}(undef, nbvec, nk_fbz)

    ikisym2ih = WannierIO.build_mapping_ik_isym(
        littlegroup_reps; nkpts_ibz = length(kpoints_ibz), n_symops = length(symops)
    )
    # Get mapping between equivalent b vectors
    b2b = get_equivalence_mappings(bvectors, symops)

    for (ikf, kf) in enumerate(kpoints_fbz)
        # isym: from ki to kf
        iki, isym_kf = fbz2ibz[ikf]
        ki = kpoints_ibz[iki]

        for (ibf, bf) in enumerate(bvectors)
            # In CPC Eq. 6, g₀(k_f) is from ki to kf, we need its inverse
            # g₀⁻¹(k_f) to go from kf to ki, to obtain b_i:
            # b_i = Wk⁻¹ b_f (with time reversal: b_i = -Wk⁻¹ b_f)
            ibi = b2b[ibf, isym_kf]
            bi = bvectors[ibi]

            # We need symmetries at ki + bi:
            # 1. find the index of ki + bi in the FBZ
            ikbi_fbz = findfirst(isequivalent(ki + bi), kpoints_fbz)
            # 2. find the index of ki+bi in the IBZ, and symmetry to get ki+bi in IBZ
            ikbi_ibz, isym_kbi = fbz2ibz[ikbi_fbz]
            # 3. find the index of kf + bf in the FBZ
            ikbf_fbz = findfirst(isequivalent(kf + bf), kpoints_fbz)
            # 4. find the index of kf+bf in the IBZ, and symmetry to get kf+bf in IBZ
            isym_kbf = fbz2ibz[ikbf_fbz][2]

            kpb_k_ibz[ibi, iki] == ikbi_ibz ||
                error("Mismatch in k+b index at iki=$(iki) ibi=$(ibi)")

            kpb_k_fbz[ibf, ikf] = ikbf_fbz
            G = kf + bf - kpoints_fbz[ikbf_fbz]
            isapprox(G, round.(G); atol = 1.0e-8) ||
                error("Non-integer G vector at ikf=$(ikf) ibf=$(ibf), G=$G")
            kpb_G_fbz[ibf, ikf] = Vec3{Int}(round.(Int, G))

            # get equivalent operation of h = g₀⁻¹(ki+bi) ∘ g₀⁻¹(kf) ∘ g₀(kf+bf)
            isym_h, factor, T = compose_symops(
                spinors, symops, [isym_kbi, isym_kf, isym_kbf], [true, true, false]
            )

            ih = ikisym2ih[ikbi_ibz][isym_h]
            d = littlegroup_reps[ih].d
            Mf_view = view(Mf, :, :, ibf, ikf)
            M_ibz_view = view(M_ibz, :, :, ibi, iki)
            # d(ĥ, k_b) passes through the antiunitary g₀(ki+bi): conjugate
            if symops[isym_kbi].time_reversal
                Mf_view .= M_ibz_view * conj.(d) * factor
            else
                Mf_view .= M_ibz_view * d * factor
            end
            # ⟨g₀(k_f)ψ| ...: overall conjugation for antiunitary g₀(k_f)
            if symops[isym_kf].time_reversal
                Mf_view .= conj.(Mf_view)
            end

            # Phase e^{-i b_f·v₀} where v₀ is the translation of g₀(k_f).
            # Identical for unitary and antiunitary g₀(k_f): the conjugation
            # from ⟨g₀ψ| combines with b_i = ∓Wk⁻¹b_f to give the same sign.
            θ1 = dot(bf, symops[isym_kf].v)
            # Phase e^{-i k_b·T} from the lattice translation t_T; it is
            # conjugated once for each antiunitary operation it passes through
            # (g₀(ki+bi), ĥ, and the final g₀(k_f) conjugation), giving a net
            # sign (-1)^(t_rev_kbi ⊻ t_rev_h ⊻ t_rev_kf) = (-1)^(t_rev_kbf).
            θ2 = dot(kpoints_ibz[ikbi_ibz], T)
            if symops[isym_kbf].time_reversal
                θ2 = -θ2
            end
            phase = exp(-im * 2π * (θ1 + θ2))
            Mf_view .*= phase
        end
    end
    iszero(kpb_k_fbz) && error("Some k+b points are not assigned in FBZ")

    return Mf, kpb_k_fbz, kpb_G_fbz
end

function reconstruct_overlaps(
        M_ibz::AbstractArray{<:Complex, 4},
        kstencil_ibz::KSpaceStencil,
        kstencil_fbz::KSpaceStencil,
        fbz2ibz::AbstractVector,
        spinors::Bool,
        symops::AbstractVector{SymOp},
        littlegroup_reps::AbstractVector{<:LittleGroupRep},
    )
    return reconstruct_overlaps(
        M_ibz,
        kstencil_ibz.kpb_k,
        kstencil_ibz.kpb_G,
        kstencil_ibz.kpoints,
        get_bvectors(kstencil_fbz; fractional = true),
        kstencil_fbz.kpoints,
        fbz2ibz,
        spinors,
        symops,
        littlegroup_reps,
    )
end

"""
    $(SIGNATURES)

Reorder the overlap matrices according to the b vector ordering in the stencil.

# Arguments
- `M`: overlap matrices at each kpoint.
- `kpb_k`: k+b index mapping of `M`.
- `kpb_G`: G-vector part of k+b of `M`.
- `kstencil`: k-space stencil defining the desired b vector ordering and kpoints.

# Return
- `M`: reordered overlap matrices at each kpoint.
"""
function reorder(
        M::AbstractArray{<:Complex, 4},
        kpb_k::AbstractMatrix{Int},
        kpb_G::AbstractMatrix,
        kstencil::KSpaceStencil,
    )
    nbvec = n_bvectors(kstencil)
    nkpts = n_kpoints(kstencil)
    nkpts == size(kpb_k, 2) || error("Mismatch in number of kpoints")
    nbvec == size(kpb_k, 1) || error("Mismatch in number of b vectors")

    M_new = zeros_overlap(eltype(M), nkpts, nbvec, size(M, 1))

    for ik in 1:nkpts
        bvecs = get_bvectors(kstencil, ik; fractional = true)
        for (ib, b) in enumerate(bvecs)
            ib0 = index_bvector(kstencil.kpoints, kpb_k, kpb_G, ik, b)
            isnothing(ib0) && error("No matching bvector found for ik=$ik ib=$ib")
            view(M_new, :, :, ib, ik) .= view(M, :, :, ib0, ik)
        end
    end

    return M_new
end
