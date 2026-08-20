using LinearAlgebra
using SparseArrays: SparseMatrixCSC, sparse, droptol!, rowvals, nonzeros, nzrange
using WannierIO: SymOp, OrbitalRep, LittleGroupRep

export SymmetryConstraint, symmetry_constraint

# -----------------------------------------------------------------------------
# Symmetry-constrained (SAWF) localization on the irreducible Brillouin zone.
#
# Implements the design of `unfold/ibz_variational_sawf.tex`: the gauge is
# optimized only at IBZ kpoints, subject to the little-group covariance
# constraint (enforced by the projector 𝒫 = `project_covariant!`), and the
# spread functional is evaluated either by
#   Level 1: expanding the gauge to the full mesh and running the standard
#            full-mesh kernels, then pulling the gradient back to the IBZ, or
#   Level 2: keeping every band-dimension product on the IBZ and reaching the
#            star members through precomputed nw×nw orbital transports
#            (`Rmat`, `Lmat`) and phases — the "transport theorem".
#
# Conventions follow src/symmetry.jl (standard Seitz; `unfold_gauge`,
# `unfold_overlaps`, `merge_symops`). Time-reversal (antiunitary) operations
# are supported through explicit conjugation flags; 𝒦_x[A] below denotes
# conj(A) if flag x is set.
#
# Key identities in code objects (A(g,k) is the nw×nw expansion matrix of
# `unfold_gauge`, i.e. Uf = 𝒦_f[Ui · A(g₀(kf), ki)]):
#   M^{(kf,bf)}  = phase · 𝒦_f[ M^{(ki,bi)} · 𝒦_bi[d(ĥ,kb)] ]     (F3)
#   MU^{(kf,bf)} = phase · 𝒦_f[ MU_i · R ],   MU_i = M^{(ki,bi)} U(ki+bi)
#   M̃^{(kf,bf)} = phase · 𝒦_f[ L† · M̃_i · R ],  M̃_i = U(ki)† MU_i
# with L = A(g₀(kf), ki) and
#   R = 𝒦_bi[ A(g₀(ki+bi), kb)† · 𝒦_h[ A(ĥ, kb)† · A(g₀(kf+bf), kb) ] ].
# The R identity uses the covariance constraint
#   U(ki) = d(ĥ,ki) · 𝒦_h[ U(ki) · A(ĥ, ki) ]   ∀ ĥ ∈ G_{ki},
# whose group average is exactly `project_covariant!`.
# -----------------------------------------------------------------------------

"""
Precomputed symmetry tables for IBZ-constrained localization.

Built once by [`symmetry_constraint`](@ref); consumed by the Level-1/Level-2
objective evaluations and by [`project_covariant!`](@ref). All matrices are
`n_wann × n_wann` orbital-space transports; band-dimension objects (the
little-group `d` matrices) are stored per (IBZ kpoint, little-group element)
for the projector, and per (FBZ kpoint, b-vector) for the Level-2 seeds.
"""
struct SymmetryConstraint{T <: Real}
    nk_fbz::Int
    nk_ibz::Int
    nbvecs::Int
    nwann::Int
    nbands::Int

    # ikf -> (iki, isym) with kf = g₀(kf) ki, g₀(kf) = symops[isym]
    fbz2ibz::Vector{Tuple{Int, Int}}
    # iki -> ikf with kpoints_fbz[ikf] == kpoints_ibz[iki] (identity map member)
    ibz2fbz::Vector{Int}
    # iki -> all star members ikf
    stars::Vector{Vector{Int}}

    # Projector data, per IBZ kpoint: (d̂(ĥ), A(ĥ, ki), trev(ĥ)) for ĥ ∈ G_ki,
    # with d̂ the little-group matrix masked on symmetry-broken bands. Stored
    # sparse: d̂ is block-diagonal over degenerate multiplets and A is block
    # monomial, so the group average costs O(nnz) instead of dense GEMMs.
    proj::Vector{Vector{Tuple{SparseMatrixCSC{Complex{T}, Int}, SparseMatrixCSC{Complex{T}, Int}, Bool}}}
    # Bands whose symmetry partners are fully inside the window (per IBZ kpoint).
    # A covariant gauge cannot carry weight on the others (their little-group
    # rows/columns are truncated), so the projector masks them out.
    band_ok::Vector{BitVector}
    # Number of group-average iterations the projector applies. Fixed at build
    # time (probed to stagnation on a test vector) so that 𝒫 is *exactly*
    # linear — an input-dependent iteration count would make the objective
    # composition Ω∘expand∘𝒫 only approximately differentiable.
    proj_niter::Int

    # Gauge expansion (C2), per FBZ kpoint: U(kf) = 𝒦_f[U(ki) · Lmat[ikf]].
    # Block monomial (permutation of site blocks × small dense blocks), hence
    # stored sparse; same for the Level-2 transports Rmat below.
    Lmat::Vector{SparseMatrixCSC{Complex{T}, Int}}
    trev_f::BitVector

    # Pass-0 tables, per (ibi, iki): U(ki+bi) = 𝒦_ib[U(kb) · Aib], kb IBZ point
    ikb::Matrix{Int}
    Aib::Matrix{Matrix{Complex{T}}}
    trev_ib::BitMatrix

    # Hermiticity-pair tables for the pass-0 halving of `_fg2_core!`:
    # `ikpb_fbz[ibi, iki]` is the FBZ index of ki + bi (the base point of the
    # dagger member (ki+bi, −bi) of the pair (ki, bi)), and `opp_b[ib]` the
    # index of −b[ib] in the global b shell. The induced partner map
    # P(iki, ibi) = (ikb[ibi, iki], ibi_of[opp_b[ibi], ikpb_fbz[ibi, iki]])
    # is generally *not* an involution (P² can land on a little-group-related
    # b at ki); `_fg2_core!` derives only the 2-cycle pairs of P and needs
    # the dagger member to star-map to the stored kb (asserted at build time).
    ikpb_fbz::Matrix{Int}
    opp_b::Vector{Int}

    # Level-2 star tables, per (ibf, ikf): b index at the IBZ source, the
    # unfolding phase (± sign of merge_symops folded in), the orbital
    # transport R, and the band-space d(ĥ, kb) with its conjugation flag
    # (needed to reconstruct full M^{(kf,bf)}; Level 2 itself only needs R)
    ibi_of::Matrix{Int}
    phase::Matrix{Complex{T}}
    Rmat::Matrix{SparseMatrixCSC{Complex{T}, Int}}
    dmat::Matrix{Matrix{Complex{T}}}
    trev_bi::BitMatrix

    # Global b shell (identical ordering at every kpoint; checked at build)
    bvec_cart::Vector{Vec3{T}}
    bweights::Vector{T}
end

n_kpoints_ibz(sc::SymmetryConstraint) = sc.nk_ibz

"""
    $(SIGNATURES)

Rebuild a k-space stencil whose neighbor lists follow the *global* b-vector
ordering of the first kpoint, i.e. `kpb_k[ib, ik]` is the kpoint index of
`k + b[ib]` for the same `b[ib]` at every `ik`. The symmetry tables of
[`symmetry_constraint`](@ref) and the outputs of [`unfold_overlaps`](@ref)
use this ordering, so full-mesh models built from them need this stencil.
"""
function globalize_stencil(kstencil::KspaceStencil)
    kpts = kstencil.kpoints
    bvecs = get_bvectors(kstencil; fractional = true)
    nk = length(kpts)
    nb = length(bvecs)
    kpb_k = zeros(Int, nb, nk)
    kpb_G = Matrix{Vec3{Int}}(undef, nb, nk)
    for ik in 1:nk, ib in 1:nb
        kpb = kpts[ik] + bvecs[ib]
        ikpb = findfirst(isequiv(kpb), kpts)
        isnothing(ikpb) && error("k+b not on the mesh: ik=$ik ib=$ib")
        G = kpb - kpts[ikpb]
        isapprox(G, round.(G); atol = 1.0e-7) || error("non-integer G at ik=$ik ib=$ib")
        kpb_k[ib, ik] = ikpb
        kpb_G[ib, ik] = Vec3{Int}(round.(Int, G))
    end
    return KspaceStencil(kstencil.recip_lattice, kpts, kpb_k, kpb_G)
end

"""
Expansion matrix ``A(g_{isym}, k)`` such that the gauge at ``k_f = g_{isym} k``
is `Uf = 𝒦[Ui * A]` (conjugated when the operation is antiunitary): the
matrix `D(g⁻¹) .* transpose(phases)` of `unfold_gauge`, with the
stored-inverse lattice mismatch folded into the translations.
"""
function _expansion_matrix(
        isym::Integer,
        k::AbstractVector,
        symops::AbstractVector{SymOp},
        orbital_reps::AbstractVector{<:OrbitalRep},
        Rs::AbstractVector,
    )
    isinv = symops[isym].isym_inv
    L = inverse_translation_mismatch(symops, isym)
    D = orbital_reps[isinv].D
    phases = [exp(-im * 2π * dot(k, Ri - L)) for Ri in Rs[isinv]]
    return Matrix{ComplexF64}(D .* transpose(phases))
end

_kconj(A::AbstractArray, trev::Bool) = trev ? conj.(A) : A
_kconj(a::Number, trev::Bool) = trev ? conj(a) : a

"""
    $(SIGNATURES)

Build the [`SymmetryConstraint`](@ref) tables from the full-mesh k stencil,
the symmetry data of an `.isym` file, and the WF centers (fractional).

# Arguments
- `kstencil`: full-mesh k-space stencil (defines kpoints and the b shell).
- `isym`: named tuple from `WannierIO.read_isym` (with `littlegroup_reps`
  rescaled by [`rescale!`](@ref) beforehand).
- `centers`: WF centers in fractional coordinates (from the `.nnkp`
  projections), used for the orbital translation vectors.
"""
function symmetry_constraint(
        kstencil::KspaceStencil{T},
        isym,
        centers::AbstractVector;
        eig_ibz::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
        atol_degeneracy::Real = 1.0e-4,
    ) where {T}
    symops = isym.symops
    kpts_ibz = isym.kpoints_ibz
    kpts_fbz = kstencil.kpoints
    orbital_reps = isym.orbital_reps
    littlegroup_reps = isym.littlegroup_reps
    spinors = isym.spinors

    nk_fbz = length(kpts_fbz)
    nk_ibz = length(kpts_ibz)
    nwann = size(orbital_reps[1].D, 1)
    nbands = size(littlegroup_reps[1].d, 1)

    Rs = find_wf_symmetry_translations(centers, symops, orbital_reps)
    f2i_raw = get_kpoint_mappings(kpts_fbz, kpts_ibz, symops)
    fbz2ibz = [Tuple(x) for x in f2i_raw]

    # Global b shell, taken from the first kpoint. All symmetry tables (and
    # the IBZ overlaps, cf. `unfold_overlaps`) use this ordering at *every*
    # kpoint; the per-kpoint orderings of a w90 stencil differ, so full-mesh
    # models consumed together with these tables must be built on
    # [`globalize_stencil`](@ref).
    bvecs_frac = get_bvectors(kstencil; fractional = true)
    nbvecs = length(bvecs_frac)
    recip = reciprocal_lattice(kstencil)
    bvec_cart = [Vec3{T}(recip * b) for b in bvecs_frac]
    bweights = Vector{T}(kstencil.bweights)

    # stars and IBZ representatives inside the FBZ mesh
    stars = [Int[] for _ in 1:nk_ibz]
    for (ikf, (iki, _)) in enumerate(fbz2ibz)
        push!(stars[iki], ikf)
    end
    ibz2fbz = zeros(Int, nk_ibz)
    for (iki, ki) in enumerate(kpts_ibz)
        ikf = findfirst(isequiv(ki), kpts_fbz)
        isnothing(ikf) && error("IBZ kpoint $iki not found in the FBZ mesh")
        ibz2fbz[iki] = ikf
    end

    # little-group (projector) tables
    ikisym2ih = WannierIO.build_mapping_ik_isym(
        littlegroup_reps; nkpts_ibz = nk_ibz, n_symops = length(symops)
    )
    CT = Complex{T}
    SM = SparseMatrixCSC{CT, Int}
    proj = Vector{Vector{Tuple{SM, SM, Bool}}}(undef, nk_ibz)
    band_ok = Vector{BitVector}(undef, nk_ibz)
    for iki in 1:nk_ibz
        entries = Tuple{Matrix{CT}, Matrix{CT}, Bool}[]
        for (isym, op) in enumerate(symops)
            ih = ikisym2ih[iki][isym]
            isnothing(ih) && continue
            d = Matrix{CT}(littlegroup_reps[ih].d)
            A = _expansion_matrix(isym, kpts_ibz[iki], symops, orbital_reps, Rs)
            push!(entries, (d, A, op.time_reversal))
        end
        # Bands whose multiplet is cut by the window have non-unit row/column
        # norms in some d(ĥ); a covariant gauge has zero weight there. Masking
        # those rows/columns makes the remaining d's exactly unitary (energy
        # blocks decouple), so one group average is idempotent.
        ok = trues(nbands)
        for (d, _, _) in entries, n in 1:nbands
            if abs(1 - sum(abs2, view(d, :, n))) > 1.0e-4 ||
                    abs(1 - sum(abs2, view(d, n, :))) > 1.0e-4
                ok[n] = false
            end
        end
        # With eigenvalues available, extend the mask to whole energy
        # multiplets: symmetry can only be broken multiplet-wise, so a
        # partially flagged degenerate cluster is masked in full (robust
        # against representation-matrix noise around the detection threshold).
        if eig_ibz !== nothing
            E = view(eig_ibz, :, iki)
            lo = 1
            for n in 1:nbands
                if n == nbands || E[n + 1] - E[n] > atol_degeneracy
                    any(!, view(ok, lo:n)) && (ok[lo:n] .= false)
                    lo = n + 1
                end
            end
        end
        proj[iki] = map(entries) do (d, A, trev)
            dm = copy(d)
            dm[.!ok, :] .= 0
            dm[:, .!ok] .= 0
            # entries below the data noise floor of the representation
            # matrices are dropped; they only add cost, not accuracy
            (droptol!(sparse(dm), 1.0e-8), droptol!(sparse(A), 1.0e-8), trev)
        end
        band_ok[iki] = ok
    end

    # Probe how many group averages reach the (data-limited) fixed point.
    probe = ComplexF64.(reshape(1:(nbands * nwann * nk_ibz), nbands, nwann, nk_ibz))
    probe ./= norm(probe)
    for iki in 1:nk_ibz
        view(probe, .!band_ok[iki], :, iki) .= 0
    end
    proj_niter = 1
    resid = Inf
    for _ in 1:30
        prev = copy(probe)
        _covariant_average_entries!(probe, proj, nk_ibz)
        r = norm(probe - prev)
        (r <= 1.0e-11 || r > 0.5 * resid) && break
        resid = r
        proj_niter += 1
    end

    # gauge expansion tables (C2)
    Lmat = Vector{SM}(undef, nk_fbz)
    trev_f = falses(nk_fbz)
    for ikf in 1:nk_fbz
        iki, isym = fbz2ibz[ikf]
        Lmat[ikf] = droptol!(
            sparse(_expansion_matrix(isym, kpts_ibz[iki], symops, orbital_reps, Rs)),
            1.0e-12,
        )
        trev_f[ikf] = symops[isym].time_reversal
    end

    # per-(ibf, ikf) star tables and per-(ibi, iki) pass-0 tables
    b2b = get_equivalence_mappings(bvecs_frac, symops)
    ikb = zeros(Int, nbvecs, nk_ibz)
    Aib = Matrix{Matrix{CT}}(undef, nbvecs, nk_ibz)
    trev_ib = falses(nbvecs, nk_ibz)
    ikpb_fbz = zeros(Int, nbvecs, nk_ibz)
    opp_b = map(bvecs_frac) do b
        ib2 = findfirst(isequiv(-b), bvecs_frac)
        isnothing(ib2) && error("b shell is not inversion-closed: -b missing for b = $b")
        ib2
    end
    opp_b[opp_b] == 1:nbvecs || error("opp_b is not an involution")
    ibi_of = zeros(Int, nbvecs, nk_fbz)
    phase = zeros(CT, nbvecs, nk_fbz)
    Rmat = Matrix{SM}(undef, nbvecs, nk_fbz)
    dmat = Matrix{Matrix{CT}}(undef, nbvecs, nk_fbz)
    trev_bi = falses(nbvecs, nk_fbz)

    for (ikf, kf) in enumerate(kpts_fbz)
        iki, isym_kf = fbz2ibz[ikf]
        ki = kpts_ibz[iki]
        for (ibf, bf) in enumerate(bvecs_frac)
            ibi = b2b[ibf, isym_kf]
            bi = bvecs_frac[ibi]

            ikbi_fbz = findfirst(isequiv(ki + bi), kpts_fbz)
            ikbi_ibz, isym_kbi = fbz2ibz[ikbi_fbz]
            kb = kpts_ibz[ikbi_ibz]
            ikbf_fbz = findfirst(isequiv(kf + bf), kpts_fbz)
            isym_kbf = fbz2ibz[ikbf_fbz][2]

            # ĥ = g₀⁻¹(ki+bi) ∘ g₀⁻¹(kf) ∘ g₀(kf+bf) ∈ G_kb  (up to t_T)
            isym_h, factor, T_lat = merge_symops(
                spinors, symops, [isym_kbi, isym_kf, isym_kbf], [true, true, false]
            )
            ih = ikisym2ih[ikbi_ibz][isym_h]
            isnothing(ih) && error("ĥ is not in the little group of kb")

            # phases exactly as in `unfold_overlaps`
            θ1 = dot(bf, symops[isym_kf].v)
            θ2 = dot(kb, T_lat)
            if symops[isym_kbf].time_reversal
                θ2 = -θ2
            end

            # orbital transport R = 𝒦_bi[ Aib† · 𝒦_h[ A_h† · A_bf ] ]
            A_h = _expansion_matrix(isym_h, kb, symops, orbital_reps, Rs)
            A_bi = _expansion_matrix(isym_kbi, kb, symops, orbital_reps, Rs)
            A_bf = _expansion_matrix(isym_kbf, kb, symops, orbital_reps, Rs)
            inner = _kconj(A_h' * A_bf, symops[isym_h].time_reversal)
            R = _kconj(A_bi' * inner, symops[isym_kbi].time_reversal)

            ibi_of[ibf, ikf] = ibi
            phase[ibf, ikf] = factor * exp(-im * 2π * (θ1 + θ2))
            Rmat[ibf, ikf] = droptol!(sparse(R), 1.0e-12)
            dmat[ibf, ikf] = _kconj(
                Matrix{CT}(littlegroup_reps[ih].d), symops[isym_kbi].time_reversal
            ) .* factor
            trev_bi[ibf, ikf] = symops[isym_kbi].time_reversal

            # pass-0 tables (fill once per (ibi, iki); star members agree)
            if ikb[ibi, iki] == 0
                ikb[ibi, iki] = ikbi_ibz
                Aib[ibi, iki] = A_bi
                trev_ib[ibi, iki] = symops[isym_kbi].time_reversal
                ikpb_fbz[ibi, iki] = ikbi_fbz
            end
        end
    end
    any(iszero, ikb) && error("Some (ibi, iki) pass-0 entries were never filled")
    # Hermiticity-pair consistency: the dagger member (ki+bi, −bi) of the
    # pair (ki, bi) must star-map to the same IBZ point kb as the pass-0
    # tables (so `_fg2_core!` can read the partner M̃ at (kb, ibi2)).
    for iki in 1:nk_ibz, ibi in 1:nbvecs
        fbz2ibz[ikpb_fbz[ibi, iki]][1] == ikb[ibi, iki] ||
            error("Hermiticity-pair map inconsistent at (ibi = $ibi, iki = $iki)")
    end

    return SymmetryConstraint{T}(
        nk_fbz, nk_ibz, nbvecs, nwann, nbands,
        fbz2ibz, ibz2fbz, stars,
        proj, band_ok, proj_niter,
        Lmat, trev_f,
        ikb, Aib, trev_ib,
        ikpb_fbz, opp_b,
        ibi_of, phase, Rmat, dmat, trev_bi,
        bvec_cart, bweights,
    )
end

# -----------------------------------------------------------------------------
# Covariance projector 𝒫 and gauge expansion / pullback
# -----------------------------------------------------------------------------

"""
    $(SIGNATURES)

One group average of the little-group action,
`𝒜[U](ki) = (1/N_h) Σ_ĥ d(ĥ) 𝒦_h[U(ki) A(ĥ, ki)]`. Self-adjoint w.r.t.
the real inner product `Re tr(X†Y)` (because `d(ĥ)† = d(ĥ⁻¹)` holds even for
truncated band windows), but idempotent only when the band window is closed
under the little group at every IBZ kpoint. When the window cuts through a
degenerate multiplet the `d` matrices are contractions on the broken rows,
and `𝒜` is a self-adjoint contraction whose *iterates* converge to the
orthogonal projector onto the covariant gauges — see
[`project_covariant!`](@ref).
"""
function _covariant_average_entries!(U_ibz, proj, nk_ibz)
    nb, nw = size(U_ibz, 1), size(U_ibz, 2)
    Uk = similar(U_ibz, nb, nw)
    UA = similar(U_ibz, nb, nw)
    acc = similar(U_ibz, nb, nw)
    for iki in 1:nk_ibz
        Uk .= view(U_ibz, :, :, iki)
        fill!(acc, 0)
        for (d, A, trev) in proj[iki]
            mul!(UA, Uk, A)
            trev && (UA .= conj.(UA))
            mul!(acc, d, UA, true, true)
        end
        view(U_ibz, :, :, iki) .= acc ./ length(proj[iki])
    end
    return U_ibz
end

covariant_average!(U_ibz::AbstractArray{<:Complex, 3}, sc::SymmetryConstraint) =
    _covariant_average_entries!(U_ibz, sc.proj, sc.nk_ibz)

"""
    $(SIGNATURES)

Covariance projector ``𝒫``: mask the symmetry-broken bands, then apply
[`covariant_average!`](@ref) a *fixed* number of times (`sc.proj_niter`,
probed at build time). Idempotent (to the numerical quality of the `d`
matrices) and self-adjoint; its image is the space of covariant gauges. The
fixed iteration count keeps 𝒫 exactly linear, so objective compositions
through 𝒫 are exactly differentiable. Overwrites `U_ibz`.
"""
function project_covariant!(
        U_ibz::AbstractArray{<:Complex, 3}, sc::SymmetryConstraint;
        niter::Integer = sc.proj_niter,
    )
    for iki in 1:sc.nk_ibz
        view(U_ibz, .!sc.band_ok[iki], :, iki) .= 0
    end
    for _ in 1:niter
        covariant_average!(U_ibz, sc)
    end
    return U_ibz
end

project_covariant(U_ibz::AbstractArray{<:Complex, 3}, sc::SymmetryConstraint; kwargs...) =
    project_covariant!(copy(U_ibz), sc; kwargs...)

"""
    $(SIGNATURES)

Residual `max_k ‖U(ki) − 𝒫[U](ki)‖` measuring how far the gauge is from
covariance (0 for a symmetry-adapted gauge).
"""
function covariance_residual(U_ibz::AbstractArray{<:Complex, 3}, sc::SymmetryConstraint)
    P = project_covariant(U_ibz, sc)
    return maximum(norm(view(U_ibz, :, :, ik) - view(P, :, :, ik)) for ik in 1:sc.nk_ibz)
end

"""
    $(SIGNATURES)

Expand IBZ gauges to the full mesh, `U(kf) = 𝒦_f[U(ki) Lmat(kf)]` (C2).
Writes into `U_fbz` and returns it.
"""
function expand_gauges!(
        U_fbz::AbstractArray{<:Complex, 3},
        U_ibz::AbstractArray{<:Complex, 3},
        sc::SymmetryConstraint,
    )
    for ikf in 1:sc.nk_fbz
        iki = sc.fbz2ibz[ikf][1]
        Uf = view(U_fbz, :, :, ikf)
        mul!(Uf, view(U_ibz, :, :, iki), sc.Lmat[ikf])
        sc.trev_f[ikf] && (Uf .= conj.(Uf))
    end
    return U_fbz
end

expand_gauges(U_ibz::AbstractArray{<:Complex, 3}, sc::SymmetryConstraint) =
    expand_gauges!(
        similar(U_ibz, size(U_ibz, 1), size(U_ibz, 2), sc.nk_fbz), U_ibz, sc
    )

"""
    $(SIGNATURES)

Pull a full-mesh canonical gradient `dΩ/dU*(kf)` back to the IBZ variables:
`G(ki) = Σ_{kf ∈ star(ki)} 𝒦_f[G(kf)] Lmat(kf)†` (adjoint of the expansion).
The covariance projector is *not* applied here.
"""
function pullback_gauges!(
        G_ibz::AbstractArray{<:Complex, 3},
        G_fbz::AbstractArray{<:Complex, 3},
        sc::SymmetryConstraint,
    )
    fill!(G_ibz, 0)
    for ikf in 1:sc.nk_fbz
        iki = sc.fbz2ibz[ikf][1]
        Gf = _kconj(view(G_fbz, :, :, ikf), sc.trev_f[ikf])
        mul!(view(G_ibz, :, :, iki), Gf, sc.Lmat[ikf]', true, true)
    end
    return G_ibz
end

"""
    $(SIGNATURES)

Extract the IBZ slices of a full-mesh gauge array (at the representative
FBZ indices of the IBZ kpoints).
"""
function extract_ibz_gauges(U_fbz::AbstractArray{<:Complex, 3}, sc::SymmetryConstraint)
    return U_fbz[:, :, sc.ibz2fbz]
end

"""
    $(SIGNATURES)

Reconstruct the full-mesh overlaps from the IBZ overlaps using the cached
tables, `M^{(kf,bf)} = phase · 𝒦_f[M^{(ki,bi)} · d]` — the same result as
[`unfold_overlaps`](@ref) but table-driven; used for validation.
"""
function unfold_overlaps_cached(
        M_ibz::AbstractArray{<:Complex, 4}, sc::SymmetryConstraint
    )
    nb = size(M_ibz, 1)
    Mf = zeros_overlap(eltype(M_ibz), sc.nk_fbz, sc.nbvecs, nb)
    for ikf in 1:sc.nk_fbz
        iki = sc.fbz2ibz[ikf][1]
        for ibf in 1:sc.nbvecs
            ibi = sc.ibi_of[ibf, ikf]
            Mv = view(Mf, :, :, ibf, ikf)
            mul!(Mv, view(M_ibz, :, :, ibi, iki), sc.dmat[ibf, ikf])
            sc.trev_f[ikf] && (Mv .= conj.(Mv))
            Mv .*= sc.phase[ibf, ikf]
        end
    end
    return Mf
end

# -----------------------------------------------------------------------------
# Level-1 evaluation: expand the IBZ gauge to the full mesh, run the standard
# full-mesh spread/gradient kernels, and pull the gradient back to the IBZ.
# -----------------------------------------------------------------------------

"""
Scratch buffers for symmetry-constrained localization. `full` carries the
standard full-mesh [`Workspace`](@ref) (its `U`/`GU`/`MU`/`UtMU` buffers are
used by the Level-1 kernels); the `*_ibz` arrays hold the IBZ variables in
the `(X, Y)` disentanglement layout.
"""
struct SymmetricWorkspace{T}
    full::Workspace{T}
    U_ibz::Array{Complex{T}, 3}
    G_ibz::Array{Complex{T}, 3}
    X_ibz::Array{Complex{T}, 3}
    Y_ibz::Array{Complex{T}, 3}
    frozen_ibz::BitMatrix
end

function SymmetricWorkspace(model::Model, sc::SymmetryConstraint{T}) where {T}
    nb, nw = n_bands(model), n_wannier(model)
    nb == sc.nbands && nw == sc.nwann ||
        error("model and symmetry constraint sizes do not match")
    n_kpoints(model) == sc.nk_fbz ||
        error("model must live on the full mesh of the symmetry constraint")
    full = Workspace(model)
    U_ibz = zeros(Complex{T}, nb, nw, sc.nk_ibz)
    G_ibz = zeros(Complex{T}, nb, nw, sc.nk_ibz)
    X_ibz = zeros(Complex{T}, nw, nw, sc.nk_ibz)
    Y_ibz = zeros(Complex{T}, nb, nw, sc.nk_ibz)
    frozen_ibz = model.frozen_bands[:, sc.ibz2fbz]
    return SymmetricWorkspace(full, U_ibz, G_ibz, X_ibz, Y_ibz, frozen_ibz)
end

"""
    $(SIGNATURES)

Level-1 fused value/gradient of the symmetry-constrained spread.

`xy` packs the `(X, Y)` blocks at the IBZ kpoints. The gauge is decoded,
projected onto the covariant subspace, expanded to the full mesh (C2), and the
standard full-mesh kernels evaluate Ω and `dΩ/dU*`; the gradient is pulled
back through the (linear, self-adjoint) expansion and projector, then packed
into `G`. The `model` must be a full-mesh model in the *global* b ordering
(see [`globalize_stencil`](@ref)) whose overlaps were unfolded from the IBZ.
"""
function symmetric_fg1!(
        F, G, xy::AbstractMatrix,
        model::Model, sc::SymmetryConstraint, ws::SymmetricWorkspace,
    )
    # decode (X,Y) -> covariant U at IBZ
    XY_to_X_Y!(ws.X_ibz, ws.Y_ibz, xy)
    X_Y_to_U!(ws.U_ibz, ws.X_ibz, ws.Y_ibz)
    project_covariant!(ws.U_ibz, sc)

    Ω = _fg1_core!(F, G === nothing ? nothing : ws.G_ibz, model, sc, ws)

    if G !== nothing
        project_covariant!(ws.G_ibz, sc)
        encode_gradient_xy!(G, ws.G_ibz, ws.X_ibz, ws.Y_ibz, ws.frozen_ibz)
    end
    return Ω
end

"""
Level-1 core: value and (unprojected) canonical gradient `dΩ/dU*(ki)` for the
covariant gauge already stored in `ws.U_ibz` — expand to the full mesh, run
the standard full-mesh kernels, pull the gradient back to the IBZ. Writes the
gradient into `G_ibz` when given. Counterpart of [`_fg2_core!`](@ref).
"""
_fg1_core!(F, G_ibz, model::Model, sc::SymmetryConstraint, ws::SymmetricWorkspace) =
    _fg1_core!(F, G_ibz, (r, _) -> r, model, sc, ws)

function _fg1_core!(
        F, G_ibz, penalty::Function,
        model::Model, sc::SymmetryConstraint, ws::SymmetricWorkspace,
    )
    expand_gauges!(ws.full.U, ws.U_ibz, sc)

    kstencil, overlaps = model.kstencil, model.overlaps
    compute_MU_UtMU!(ws.full, kstencil, overlaps, ws.full.U)

    if G_ibz !== nothing
        omega_grad!(penalty, ws.full.GU, ws.full, kstencil, overlaps)
        pullback_gauges!(G_ibz, ws.full.GU, sc)
    end

    F === nothing && return nothing
    # value of the *plain* spread; penalty terms are the objective's to add
    # (the WF centers are left in `ws.full.r` for it)
    return omega!(ws.full, kstencil, overlaps).Ω
end

# -----------------------------------------------------------------------------
# Level-2 evaluation: all band-dimension products stay on the IBZ; star
# members are reached through the nw×nw orbital transports (transport theorem)
#   M̃^{(kf,bf)} = phase · 𝒦_f[ L† M̃_i R ],
# and the gradient seeds are accumulated back onto the IBZ pairs,
#   𝒦^{(ki,bi)} = Σ_{kf ∈ star} (4 w_b/N_k) 𝒦_f[phase] · R 𝒦_f[diag T] L†,
#   dΩ/dU*(ki) = Σ_bi MU_i 𝒦^{(ki,bi)}  (then projected).
# -----------------------------------------------------------------------------

"""
Scratch buffers for the Level-2 evaluation. All arrays are IBZ-sized except
the transported diagonals `tdiag`, which live on the full mesh but hold only
`n_wann` numbers per (kpoint, b-vector).
"""
struct SymmetricWorkspace2{T}
    U_ibz::Array{Complex{T}, 3}
    G_ibz::Array{Complex{T}, 3}
    X_ibz::Array{Complex{T}, 3}
    Y_ibz::Array{Complex{T}, 3}
    frozen_ibz::BitMatrix
    # U(ki+bi), MU_i = M^{(ki,bi)} U(ki+bi), M̃_i = U(ki)† MU_i, per (ibi, iki)
    Ukb::Array{Complex{T}, 4}
    MU::Array{Complex{T}, 4}
    Mt::Array{Complex{T}, 4}
    # per-(ibi, iki) seed accumulators 𝒦
    K::Array{Complex{T}, 4}
    # transported diagonals t_n^{(kf, bf)} on the full mesh
    tdiag::Array{Complex{T}, 3}
    r::Vector{Vec3{T}}
    # fixed guiding centers for Im-log branches (zeros = principal branch)
    guide::Vector{Vec3{T}}
    # nw×nw scratch
    tmp1::Matrix{Complex{T}}
    tmp2::Matrix{Complex{T}}
end

function SymmetricWorkspace2(
        eig::AbstractMatrix, frozen_bands::AbstractMatrix{Bool}, sc::SymmetryConstraint{T}
    ) where {T}
    nb, nw, nki, nkf, nbv = sc.nbands, sc.nwann, sc.nk_ibz, sc.nk_fbz, sc.nbvecs
    size(eig, 1) == nb || error("eig has wrong number of bands")
    CT = Complex{T}
    return SymmetricWorkspace2(
        zeros(CT, nb, nw, nki), zeros(CT, nb, nw, nki),
        zeros(CT, nw, nw, nki), zeros(CT, nb, nw, nki),
        BitMatrix(frozen_bands[:, sc.ibz2fbz]),
        zeros(CT, nb, nw, nbv, nki), zeros(CT, nb, nw, nbv, nki),
        zeros(CT, nw, nw, nbv, nki), zeros(CT, nw, nw, nbv, nki),
        zeros(CT, nw, nbv, nkf), zeros(Vec3{T}, nw), zeros(Vec3{T}, nw),
        zeros(CT, nw, nw), zeros(CT, nw, nw),
    )
end

"""
    $(SIGNATURES)

Level-2 fused value/gradient of the symmetry-constrained spread, consuming
only the IBZ overlaps `M_ibz` (global b ordering, as in the `.immn` file).
Same variables and same value/gradient as [`symmetric_fg1!`](@ref), evaluated
without ever forming band-dimension objects on the full mesh.
"""
function symmetric_fg2!(
        F, G, xy::AbstractMatrix,
        M_ibz::AbstractArray{<:Complex, 4}, sc::SymmetryConstraint{T},
        ws::SymmetricWorkspace2{T},
    ) where {T}
    # decode -> covariant U at IBZ
    XY_to_X_Y!(ws.X_ibz, ws.Y_ibz, xy)
    X_Y_to_U!(ws.U_ibz, ws.X_ibz, ws.Y_ibz)
    project_covariant!(ws.U_ibz, sc)

    Ω = _fg2_core!(F, G === nothing ? nothing : ws.G_ibz, M_ibz, sc, ws)

    if G !== nothing
        project_covariant!(ws.G_ibz, sc)
        encode_gradient_xy!(G, ws.G_ibz, ws.X_ibz, ws.Y_ibz, ws.frozen_ibz)
    end
    return Ω
end

"""
Level-2 core: value and (unprojected) canonical gradient `dΩ/dU*(ki)` for the
covariant gauge already stored in `ws.U_ibz`. Writes the gradient into
`G_ibz` when given. The `penalty` function hooks a WF-center penalty into the
gradient seeds' `q_n`, exactly as in the penalty-aware `omega_grad!` of
src/spread.jl (the value of the penalty term is the objective's to add; the
WF centers are left in `ws.r` for it). The identity default is a no-op.
"""
_fg2_core!(
    F, G_ibz,
    M_ibz::AbstractArray{<:Complex, 4}, sc::SymmetryConstraint,
    ws::SymmetricWorkspace2,
) = _fg2_core!(F, G_ibz, (r, _) -> r, M_ibz, sc, ws)

function _fg2_core!(
        F, G_ibz, penalty::Function,
        M_ibz::AbstractArray{<:Complex, 4}, sc::SymmetryConstraint{T},
        ws::SymmetricWorkspace2{T},
    ) where {T}
    nw, nki, nkf, nbv = sc.nwann, sc.nk_ibz, sc.nk_fbz, sc.nbvecs
    wb = sc.bweights

    # Pass 0 (heavy, IBZ only): U(ki+bi), MU_i, M̃_i — with Hermiticity-pair
    # halving of the band-dimension GEMMs. The pair (ki, bi) has its dagger
    # member at the FBZ point ikpb = sc.ikpb_fbz[ibi, iki] with b-vector −bi:
    # in the fixed gauge M̃^{(ikpb, −bi)} = (M̃_i^{(iki,ibi)})† exactly. The
    # star tables at (ib2, ikpb), ib2 = sc.opp_b[ibi], map that same object
    # to the IBZ pair p = (kb, ibi2) = (sc.ikb[ibi,iki], sc.ibi_of[ib2,ikpb]):
    #   (M̃_q)† = phase2 · 𝒦_2[ L2† · M̃_p · R2 ],   q = (iki, ibi),
    # with phase2 = sc.phase[ib2,ikpb], L2 = sc.Lmat[ikpb],
    # R2 = sc.Rmat[ib2,ikpb], 𝒦_2 flag sc.trev_f[ikpb]. Solving (L2, R2
    # unitary) gives the partner M̃ as a cheap nw-dimensional conjugation:
    #   M̃_q = conj(phase2) · 𝒦_2[ R2† · M̃_p† · L2 ]                     (HP)
    #
    # The partner map q ↦ p is not an involution: P² can land on a
    # little-group-related b at ki. Where q and p form a true 2-cycle
    # (P(p) = q) the stored overlaps are Hermiticity-consistent and (HP) is
    # machine-exact (max 5e-12 on Si2_hse, 3e-13 on Ge4Ru4); along P-chains
    # the two data entries are related only through a little-group rotation,
    # so (HP) holds only to the data's symmetry noise (4e-7 on Ge4Ru4, 7e-4
    # on Si2_hse). Only 2-cycle pairs are therefore *derived* (the member
    # with the larger linear key (iki−1)·nbv + ibi); everything else —
    # including self-paired pairs — is *canonical* and computed directly,
    # which keeps the Level-2 value/gradient identical to Level 1 instead of
    # only data-noise-close.
    pair_key(jki, jbi) = (jki - 1) * nbv + jbi
    partner_of(jki, jbi) =
        (sc.ikb[jbi, jki], sc.ibi_of[sc.opp_b[jbi], sc.ikpb_fbz[jbi, jki]])
    function pair_derived(jki, jbi)
        kb, jbi2 = partner_of(jki, jbi)
        return pair_key(kb, jbi2) < pair_key(jki, jbi) &&
            partner_of(kb, jbi2) == (jki, jbi)
    end

    # Loop A (all pairs): U(ki+bi); and MU_i when the gradient is requested —
    # the gradient assembly needs MU for every pair, so there is no halving
    # there. The value-only path skips MU here entirely and skips U(ki+bi)
    # for the derived pairs (unused there) — that is where the pass-0 gain
    # lives.
    for iki in 1:nki, ibi in 1:nbv
        G_ibz === nothing && pair_derived(iki, ibi) && continue
        ikb = sc.ikb[ibi, iki]
        Uk = view(ws.Ukb, :, :, ibi, iki)
        mul!(Uk, view(ws.U_ibz, :, :, ikb), sc.Aib[ibi, iki])
        sc.trev_ib[ibi, iki] && (Uk .= conj.(Uk))
        if G_ibz !== nothing
            mul!(view(ws.MU, :, :, ibi, iki), view(M_ibz, :, :, ibi, iki), Uk)
        end
    end
    # Loop B (canonical pairs only): M̃_i = U(ki)† (M·U(ki+bi)); on the
    # value-only path the M·U product is formed here, for the canonical
    # pairs only, reusing the MU slice as scratch.
    for iki in 1:nki, ibi in 1:nbv
        pair_derived(iki, ibi) && continue
        MUv = view(ws.MU, :, :, ibi, iki)
        if G_ibz === nothing
            mul!(MUv, view(M_ibz, :, :, ibi, iki), view(ws.Ukb, :, :, ibi, iki))
        end
        mul!(view(ws.Mt, :, :, ibi, iki), view(ws.U_ibz, :, :, iki)', MUv)
    end
    # Loop C (derived pairs): M̃ from the canonical partner via the transport
    # (HP) — no band-dimension GEMM. The partner of a derived pair is the
    # smaller-key member of a 2-cycle, hence canonical and already computed.
    for iki in 1:nki, ibi in 1:nbv
        pair_derived(iki, ibi) || continue
        kb, ibi2 = partner_of(iki, ibi)
        ikpb = sc.ikpb_fbz[ibi, iki]
        ib2 = sc.opp_b[ibi]
        Mq = view(ws.Mt, :, :, ibi, iki)
        ws.tmp1 .= view(ws.Mt, :, :, ibi2, kb)'
        mul!(ws.tmp2, sc.Rmat[ib2, ikpb]', ws.tmp1)
        mul!(Mq, ws.tmp2, sc.Lmat[ikpb])
        sc.trev_f[ikpb] && (Mq .= conj.(Mq))
        Mq .*= conj(sc.phase[ib2, ikpb])
    end

    # Sweep 1 (light, full mesh): transported diagonals and centers.
    # `ws.guide` fixes the Im-log branches (cf. spread.jl; zeros by default).
    rg = ws.guide
    fill!(ws.r, zero(Vec3{T}))
    for ikf in 1:nkf
        iki = sc.fbz2ibz[ikf][1]
        L = sc.Lmat[ikf]
        Lrv, Lnz = rowvals(L), nonzeros(L)
        for ibf in 1:nbv
            ibi = sc.ibi_of[ibf, ikf]
            R = sc.Rmat[ibf, ikf]
            Rrv, Rnz = rowvals(R), nonzeros(R)
            Mt = view(ws.Mt, :, :, ibi, iki)
            ph = sc.phase[ibf, ikf]
            trev = sc.trev_f[ikf]
            # t_n = (L e_n)† M̃_i (R e_n): O(s²) per diagonal via the sparse
            # columns of the block-monomial transports
            @inbounds for n in 1:nw
                t = zero(Complex{T})
                for iL in nzrange(L, n), iR in nzrange(R, n)
                    t += conj(Lnz[iL]) * Mt[Lrv[iL], Rrv[iR]] * Rnz[iR]
                end
                trev && (t = conj(t))
                t *= ph
                ws.tdiag[n, ibf, ikf] = t
                gl = imaglog_guided(t, sc.bvec_cart[ibf] ⋅ rg[n])
                ws.r[n] -= gl * (wb[ibf] * sc.bvec_cart[ibf])
            end
        end
    end
    ws.r ./= nkf

    Ω = nothing
    if F !== nothing
        Ω = zero(T)
        @inbounds for ikf in 1:nkf, ibf in 1:nbv, n in 1:nw
            t = ws.tdiag[n, ibf, ikf]
            gl = imaglog_guided(t, sc.bvec_cart[ibf] ⋅ rg[n])
            Ω += wb[ibf] * (1 - abs2(t) + gl^2)
        end
        Ω = Ω / nkf - sum(r -> sum(abs2, r), ws.r)
    end

    if G_ibz !== nothing
        # Sweep 2: accumulate the seeds 𝒦 on the IBZ pairs
        fill!(ws.K, 0)
        Tn = zeros(Complex{T}, nw)
        for ikf in 1:nkf
            iki = sc.fbz2ibz[ikf][1]
            L = sc.Lmat[ikf]
            Lrv, Lnz = rowvals(L), nonzeros(L)
            trev = sc.trev_f[ikf]
            for ibf in 1:nbv
                ibi = sc.ibi_of[ibf, ikf]
                c = 4 * wb[ibf] / nkf
                # 𝒦_f[phase · diag(T)]: both factors conjugated for
                # time-reversal star members, neither otherwise
                ph = _kconj(sc.phase[ibf, ikf], trev)
                @inbounds for n in 1:nw
                    t = ws.tdiag[n, ibf, ikf]
                    q = imaglog_guided(t, sc.bvec_cart[ibf] ⋅ rg[n]) +
                        sc.bvec_cart[ibf] ⋅ penalty(ws.r[n], n)
                    s = -im * q / t - conj(t)
                    trev && (s = conj(s))
                    Tn[n] = c * ph * s
                end
                # 𝒦 += R diag(Tn) L† = Σ_n Tn[n] (R e_n)(L e_n)†:
                # O(s²) outer products over the sparse columns
                R = sc.Rmat[ibf, ikf]
                Rrv, Rnz = rowvals(R), nonzeros(R)
                Kv = view(ws.K, :, :, ibi, iki)
                @inbounds for n in 1:nw
                    for iR in nzrange(R, n)
                        rt = Rnz[iR] * Tn[n]
                        m = Rrv[iR]
                        for iL in nzrange(L, n)
                            Kv[m, Lrv[iL]] += rt * conj(Lnz[iL])
                        end
                    end
                end
            end
        end
        # assembly: dΩ/dU*(ki) = Σ_bi MU_i 𝒦
        fill!(G_ibz, 0)
        for iki in 1:nki, ibi in 1:nbv
            mul!(
                view(G_ibz, :, :, iki), view(ws.MU, :, :, ibi, iki),
                view(ws.K, :, :, ibi, iki), true, true,
            )
        end
    end

    return Ω
end

"""
    $(SIGNATURES)

Symmetry-breaking force at a covariant IBZ gauge: the relative size of the
component of the (unconstrained, IBZ-pulled-back) spread gradient that the
covariance projector removes, `‖G − 𝒫[G]‖ / ‖G‖`. Zero would mean the
constrained optimum is also an unconstrained stationary point; a sizable
value quantifies how strongly the free MLWFs want to break the imposed
symmetry (`Ω` is not invariant under the symmetry action on arbitrary
gauges, so this is generically nonzero even at convergence).
"""
function symmetry_breaking_force(
        U_ibz::AbstractArray{<:Complex, 3},
        M_ibz::AbstractArray{<:Complex, 4},
        sc::SymmetryConstraint,
        ws::SymmetricWorkspace2,
    )
    ws.U_ibz .= U_ibz
    _fg2_core!(nothing, ws.G_ibz, M_ibz, sc, ws)
    G = copy(ws.G_ibz)
    P = project_covariant(G, sc)
    return norm(G - P) / norm(G)
end

# -----------------------------------------------------------------------------
# Schur (per-irrep block) parametrization — "Phase B" of the design doc.
#
# At each IBZ kpoint the unitary little subgroup acts on the (unmasked) band
# space by d̂(ĥ) and on the orbital space by W(ĥ) = A(ĥ)†. Simultaneous
# block diagonalization (via eigenspaces of a generic commutant element)
# splits both into irrep copies; aligning every copy to a per-class reference
# irrep turns the covariance constraint into Schur block form
#   U(ki) = Σ_λ Σ_{j=1..dλ} Bb_λ[j] · C_λ · Bo_λ[j]†,
# with C_λ (m_b × m_o) free up to the Stiefel + frozen conditions, which are
# handled by the per-block DLL (X, Y) parametrization. Unitary-subgroup
# covariance is then exact by construction.
#
# Anti-unitary little-group elements â₀ (action ρ(â₀)[U] = d_a · conj(U A_a))
# permute the unitary-subgroup irrep classes, so the fixed-point condition
# U = ρ(â₀)[U] is block-wise antilinear ("corepresentation Schur blocks"):
# for every class λ there is a partner class λ′ and unitary intertwiners
# (S_b, S_o) on the copy spaces — Kronecker factors of the stacked-basis
# intertwiners, by Schur's lemma — such that
#   C_{λ′} = S_b · conj(C_λ) · S_o†.
# Pairing type (λ′ ≠ λ): C_{λ′} is derived from C_λ in the decode (S_b/S_o
# absorbed into the λ′ bases), halving those parameters. Self-paired classes
# (λ′ = λ) satisfy an antilinear block involution K[C] = S_b conj(C) S_o†
# with K² = ω = ±1 (Wigner type, S conj(S) = ω·1 with consistent signs):
#   real type (ω = +1): Takagi-factorize S = W Wᵀ and absorb W into the
#     bases — the block becomes conj(C̃) = C̃, i.e. REAL (X, Y) parameters;
#   quaternionic type (ω = −1): Youla-factorize S = W J Wᵀ (J the standard
#     interleaved symplectic form) and absorb W — the block becomes
#     conj(C̃) = J_b C̃ J_o†, i.e. quaternion-structured (X, Y) stored by
#     their independent quaternion components.
# The flat parameter vector `x` is therefore REAL: unconstrained blocks
# store re/im pairs, real-type blocks store real entries, quaternionic
# blocks store the (a, b) components of q = [a b; -conj(b) conj(a)].
# Anti-unitary covariance is exact by construction wherever the numerical
# classification succeeds; on any failure a class falls back to the soft
# 2-term coset average `_aavg` in the decode.
# -----------------------------------------------------------------------------

using Random: MersenneTwister, randn!

"""One irrep class at one IBZ kpoint: partner-stacked band/orbital bases."""
struct SchurBlock{T}
    dim::Int   # irrep dimension dλ
    mb::Int    # band multiplicity (frozen copies first)
    mo::Int    # orbital multiplicity
    mf::Int    # frozen band multiplicity
    Bb::Vector{Matrix{Complex{T}}}   # dλ matrices, n_bands × mb
    Bo::Vector{Matrix{Complex{T}}}   # dλ matrices, n_wann  × mo
    # anti-unitary handling (see `_classify_aop`): 0 = none / soft coset
    # average fallback; 1 = pairing source (free block, partner derived);
    # 2 = pairing derived (no parameters; S_b/S_o absorbed into Bb/Bo, the
    # block value is conj(C) of the source block `partner`); 3 = Wigner real
    # type (Takagi W absorbed, real block); 4 = Wigner quaternionic type
    # (Youla W absorbed, quaternion-structured block)
    akind::Int8
    partner::Int   # class index of the pairing partner (akind 1/2), else 0
end

SchurBlock{T}(dim, mb, mo, mf, Bb, Bo) where {T} =
    SchurBlock{T}(dim, mb, mo, mf, Bb, Bo, Int8(0), 0)

"""
Schur-adapted bases for all IBZ kpoints, plus the anti-unitary coset
representative (band rep, orbital matrix) where the little group is magnetic
and some class needs the soft coset-average fallback. `nx` counts REAL
parameters (see `_block_nparams`).
"""
struct SchurBasis{T}
    blocks::Vector{Vector{SchurBlock{T}}}
    aop::Vector{Union{Nothing, Tuple{Matrix{Complex{T}}, Matrix{Complex{T}}}}}
    nx::Int
end

# eigenvalue clustering: indices grouped by gaps
function _cluster(E::AbstractVector{<:Real}; rtol = 1.0e-4)
    tol = rtol * (E[end] - E[1] + 1)
    groups = Vector{UnitRange{Int}}()
    lo = 1
    for i in 1:(length(E) - 1)
        if E[i + 1] - E[i] > tol
            push!(groups, lo:i)
            lo = i + 1
        end
    end
    push!(groups, lo:length(E))
    return groups
end

# irrep copies of the unitary representation ρs (list of dim×dim matrices):
# eigenspaces of a generic commutant element. Returns (Q, ρQ) pairs.
function _irrep_copies(ρs::Vector{Matrix{CT}}, rng) where {CT}
    dim = size(ρs[1], 1)
    H = randn(rng, CT, dim, dim)
    H = Matrix(Hermitian(H + H'))
    Hc = zeros(CT, dim, dim)
    for ρ in ρs
        Hc .+= ρ * H * ρ'
    end
    Hc ./= length(ρs)
    E, V = eigen(Hermitian((Hc + Hc') / 2))
    copies = Tuple{Matrix{CT}, Vector{Matrix{CT}}}[]
    for g in _cluster(E)
        Q = V[:, g]
        ρQ = [Q' * ρ * Q for ρ in ρs]
        push!(copies, (Q, ρQ))
    end
    return copies
end

# Schur intertwiner test/alignment: S = mean(ρQ Z ρref†). Returns the unitary
# aligner (Q ← Q·u makes ρQ ≡ ρref) or nothing if inequivalent.
function _align_to(ρQ::Vector{Matrix{CT}}, ρref::Vector{Matrix{CT}}, Z) where {CT}
    d = size(Z, 1)
    (size(ρQ[1], 1) == d && size(ρref[1], 1) == d) || return nothing
    S = zeros(CT, d, d)
    for (a, b) in zip(ρQ, ρref)
        S .+= a * Z * b'
    end
    S ./= length(ρQ)
    norm(S) < 1.0e-2 * norm(Z) && return nothing
    return orthonorm_lowdin(S)
end

# Kronecker factors of an irrep-class intertwiner T ≈ kron(R, S), with R
# (d × d) acting on the partner index and S (m2 × m) on the copy index —
# exact by Schur's lemma since every copy carries the same reference irrep.
function _kron_factor(T::AbstractMatrix, d::Int, m2::Int, m::Int)
    blk(j2, j) = T[((j2 - 1) * m2 + 1):(j2 * m2), ((j - 1) * m + 1):(j * m)]
    nrm, j2m, jm = -1.0, 1, 1
    for j2 in 1:d, j in 1:d
        n = norm(blk(j2, j))
        n > nrm && ((nrm, j2m, jm) = (n, j2, j))
    end
    S = blk(j2m, jm) .* (sqrt(m) / nrm)
    R = [sum(conj.(S) .* blk(j2, j)) / m for j2 in 1:d, j in 1:d]
    return R, S, norm(kron(R, S) - T)
end

# standard interleaved symplectic form: blkdiag([0 1; -1 0], …), n even
function _jmat(::Type{CT}, n::Int) where {CT}
    J = zeros(CT, n, n)
    for i in 1:2:n
        J[i, i + 1] = 1
        J[i + 1, i] = -1
    end
    return J
end

# Takagi factor of a unitary SYMMETRIC S: unitary W with W Wᵀ = S, via the
# principal square root with the branch cut placed in the largest spectral gap
function _takagi_unitary(S::Matrix{CT}) where {CT}
    n = size(S, 1)
    n == 0 && return zeros(CT, 0, 0)
    θ = sort(angle.(eigvals(S)))
    gaps = [(i == n ? θ[1] + 2π : θ[i + 1]) - θ[i] for i in 1:n]
    i = argmax(gaps)
    α = θ[i] + gaps[i] / 2 - π
    return Matrix(sqrt(Matrix(cis(-α) .* S)) .* cis(α / 2))
end

# Youla factor of a unitary ANTISYMMETRIC S (n even): unitary W with
# W J Wᵀ = S, built two columns at a time from w₂ = -S conj(w₁)
function _youla_unitary(S::Matrix{CT}) where {CT}
    n = size(S, 1)
    W = zeros(CT, n, 0)
    for _ in 1:(n ÷ 2)
        k = argmin(vec(sum(abs2, W; dims = 2)))   # least-covered canonical axis
        v = -W * conj.(W[k, :])
        v[k] += 1
        v ./= norm(v)
        w = -S * conj.(v)
        w .-= W * (W' * w)
        w .-= v .* (v' * w)
        w ./= norm(w)
        W = hcat(W, v, w)
    end
    return W
end

# quaternion-structured (2p)×(2q) matrix M (conj(M) = J M J† in the
# interleaved convention) from its quaternion components q = [a b; -b̄ ā]
function _quat_assemble(A::AbstractMatrix{CT}, B::AbstractMatrix{CT}) where {CT}
    p, q = size(A)
    M = zeros(CT, 2p, 2q)
    for k in 1:q, i in 1:p
        M[2i - 1, 2k - 1] = A[i, k]
        M[2i - 1, 2k] = B[i, k]
        M[2i, 2k - 1] = -conj(B[i, k])
        M[2i, 2k] = conj(A[i, k])
    end
    return M
end

# inverse of `_quat_assemble` (with symmetrization = structure projection)
function _quat_extract(M::AbstractMatrix)
    p, q = size(M, 1) ÷ 2, size(M, 2) ÷ 2
    A = [(M[2i - 1, 2k - 1] + conj(M[2i, 2k])) / 2 for i in 1:p, k in 1:q]
    B = [(M[2i - 1, 2k] - conj(M[2i, 2k - 1])) / 2 for i in 1:p, k in 1:q]
    return A, B
end

# Kramers partner J·conj(v) in the interleaved convention
function _jconj(v::AbstractVector)
    w = similar(v)
    for i in 1:2:length(v)
        w[i] = conj(v[i + 1])
        w[i + 1] = -conj(v[i])
    end
    return w
end

# structured orthonormal basis (ncols columns, in Kramers pairs) of the
# column range of a quaternion-structured matrix M
function _quat_range(M::Matrix{CT}, ncols::Int) where {CT}
    n = size(M, 1)
    W = zeros(CT, n, ncols)
    nw = 0
    R = copy(M)
    while nw < ncols
        k = argmax(vec(sum(abs2, R; dims = 1)))
        v = R[:, k] ./ norm(R[:, k])
        w = -_jconj(v)   # col 2k = -J conj(col 2k-1), the `_quat_assemble` pairing
        Wv = view(W, :, 1:nw)
        w .-= Wv * (Wv' * w)
        w .-= v .* (v' * w)
        w ./= norm(w)
        W[:, nw + 1] .= v
        W[:, nw + 2] .= w
        nw += 2
        Wv = view(W, :, 1:nw)
        R .= M .- Wv * (Wv' * M)
    end
    return W
end

# Corepresentation classification of the Schur classes at one IBZ kpoint.
#
# The anti-unitary coset representative â₀ maps the isotypic component of
# class c antilinearly onto that of a partner class c′: on the stacked bases,
#   d_a · conj([Bb_c[1] … Bb_c[d]]) = [Bb_c′[1] … Bb_c′[d]] · kron(R_b, S_b),
#   A_aᵀ · conj([Bo_c[1] … Bo_c[d]]) = [Bo_c′[1] … Bo_c′[d]] · kron(R_o, S_o),
# with R_b R_o† = φ·1 (the φ is folded into S_b below). The fixed point
# U = ρ(â₀)[U] then reads block-wise C_{c′} = S_b · conj(C_c) · S_o†.
#
# Pairing type (c′ ≠ c): the lower-index class of the pair keeps its free
# parameters (akind = 1); the partner block becomes derived (akind = 2) with
# S_b/S_o absorbed into its bases, so its decode is simply
# Σ_j Bb′[j] · conj(C_c) · Bo′[j]†.
#
# Self-paired classes (c′ = c): Wigner classification by the sign ω of
# S conj(S) = ω·1. Real type (ω = +1): Takagi factors W (per frozen block on
# the band side) are absorbed into the bases, making the block constraint
# conj(C̃) = C̃ (akind = 3). Quaternionic type (ω = −1): Youla factors W with
# the interleaved symplectic J are absorbed, making the constraint
# conj(C̃) = J_b C̃ J_o† (akind = 4).
#
# Every derived identity is verified numerically (tolerance `rtol`, loose
# enough for noisy isym data); on any failure the classes fall back to
# akind = 0 (soft coset average).
function _classify_aop(
        blks::Vector{SchurBlock{T}}, da::Matrix{Complex{T}},
        Aa::Matrix{Complex{T}}, rng; rtol = 1.0e-3,
    ) where {T}
    CT = Complex{T}
    ncl = length(blks)
    SBb = [hcat(b.Bb...) for b in blks]
    SBo = [hcat(b.Bo...) for b in blks]
    Yb = [da * conj.(S) for S in SBb]
    Yo = [transpose(Aa) * conj.(S) for S in SBo]

    # intertwiner pair (S_b, S_o) mapping class c onto class c2, polished to
    # exact unitarity and to the exact frozen split of S_b; nothing if the
    # numerical Kronecker/phase structure fails
    function block_maps(c, c2)
        b, b2 = blks[c], blks[c2]
        (b.dim == b2.dim && b.mb == b2.mb && b.mo == b2.mo && b.mf == b2.mf) ||
            return nothing
        Tb = SBb[c2]' * Yb[c]
        To = SBo[c2]' * Yo[c]
        tolb, tolo = rtol * sqrt(b.dim * b.mb), rtol * sqrt(b.dim * b.mo)
        (norm(Tb' * Tb - I) < tolb && norm(To' * To - I) < tolo) || return nothing
        Rb, Sb, resb = _kron_factor(Tb, b.dim, b2.mb, b.mb)
        Ro, So, reso = _kron_factor(To, b.dim, b2.mo, b.mo)
        (resb < tolb && reso < tolo) || return nothing
        φ = tr(Rb * Ro') / b.dim
        (norm(Rb * Ro' - φ * I) < rtol * sqrt(b.dim) && abs(abs(φ) - 1) < rtol) ||
            return nothing
        # â₀ preserves the frozen (energy-window) subspace, so S_b is exactly
        # frozen-block-diagonal; enforce the split and polish per block
        mf = b.mf
        r = (mf + 1):b.mb
        (norm(Sb[1:mf, r]) + norm(Sb[r, 1:mf])) < tolb || return nothing
        Sb2 = zeros(CT, b.mb, b.mb)
        Sb2[1:mf, 1:mf] .= orthonorm_lowdin(Sb[1:mf, 1:mf])
        Sb2[r, r] .= orthonorm_lowdin(Sb[r, r])
        Sb2 .*= φ / abs(φ)
        return Sb2, Matrix{CT}(orthonorm_lowdin(So))
    end

    out = copy(blks)
    assigned = falses(ncl)
    for c in 1:ncl
        assigned[c] && continue
        assigned[c] = true
        # partner class: the antilinear image of the copy spaces of c must
        # lie in the span of exactly one class (same on band/orbital sides)
        reso = [norm(SBo[c2] * (SBo[c2]' * Yo[c]) - Yo[c]) / norm(Yo[c]) for c2 in 1:ncl]
        c2 = argmin(reso)
        resb = norm(SBb[c2] * (SBb[c2]' * Yb[c]) - Yb[c]) / norm(Yb[c])
        (reso[c2] < rtol && resb < rtol) || continue
        if c2 == c
            # Wigner classification of the self-paired antilinear involution
            fw = block_maps(c, c)
            fw === nothing && continue
            Sb, So = fw
            b = blks[c]
            mf = b.mf
            r = (mf + 1):b.mb
            tolb, tolo = rtol * sqrt(b.mb), rtol * sqrt(b.mo)
            ωb = real(tr(Sb * conj.(Sb))) / b.mb
            ωo = real(tr(So * conj.(So))) / b.mo
            (abs(abs(ωb) - 1) < rtol && abs(ωb - ωo) < rtol &&
                norm(Sb * conj.(Sb) - ωb * I) < tolb &&
                norm(So * conj.(So) - ωo * I) < tolo) || continue
            Wb = zeros(CT, b.mb, b.mb)
            local Wo, kind, Ct
            if ωb > 0   # real type: S = W Wᵀ
                Wb[1:mf, 1:mf] .= _takagi_unitary(Sb[1:mf, 1:mf])
                Wb[r, r] .= _takagi_unitary(Sb[r, r])
                Wo = _takagi_unitary(So)
                (norm(Wb * transpose(Wb) - Sb) < tolb &&
                    norm(Wo * transpose(Wo) - So) < tolo) || continue
                kind = Int8(3)
                Ct = Matrix{CT}(randn(rng, T, b.mb, b.mo))
            else        # quaternionic type: S = W J Wᵀ
                (iseven(mf) && iseven(b.mb) && iseven(b.mo)) || continue
                Wb[1:mf, 1:mf] .= _youla_unitary(Sb[1:mf, 1:mf])
                Wb[r, r] .= _youla_unitary(Sb[r, r])
                Wo = _youla_unitary(So)
                (norm(Wb * _jmat(CT, b.mb) * transpose(Wb) - Sb) < tolb &&
                    norm(Wo * _jmat(CT, b.mo) * transpose(Wo) - So) < tolo) ||
                    continue
                kind = Int8(4)
                Ct = _quat_assemble(
                    randn(rng, CT, b.mb ÷ 2, b.mo ÷ 2),
                    randn(rng, CT, b.mb ÷ 2, b.mo ÷ 2),
                )
            end
            Bb2 = [B * Wb for B in b.Bb]
            Bo2 = [B * Wo for B in b.Bo]
            # end-to-end: a random structured block value must give an
            # exactly â₀-invariant contribution
            D = sum(Bb2[j] * Ct * Bo2[j]' for j in 1:b.dim)
            norm(da * conj.(D * Aa) - D) < rtol * norm(D) || continue
            out[c] = SchurBlock{T}(b.dim, b.mb, b.mo, b.mf, Bb2, Bo2, kind, 0)
            continue
        end
        assigned[c2] && continue
        fw = block_maps(c, c2)
        bw = block_maps(c2, c)
        assigned[c2] = true
        # the two directions must invert each other (ρ(â₀)² acts trivially
        # on covariant gauges) up to a common phase μ on both sides — each
        # (S_b, S_o) pair is itself only defined up to a common phase
        fw !== nothing && bw !== nothing || continue
        Pb = bw[1] * conj.(fw[1])
        Po = bw[2] * conj.(fw[2])
        μ = tr(Pb) / blks[c].mb
        (norm(Pb - μ * I) < rtol * sqrt(blks[c].mb) &&
            norm(Po - μ * I) < rtol * sqrt(blks[c].mo) &&
            abs(abs(μ) - 1) < rtol) || continue
        b, b2 = blks[c], blks[c2]
        Bb2 = [B * fw[1] for B in b2.Bb]
        Bo2 = [B * fw[2] for B in b2.Bo]
        # end-to-end check on a random block value C: the â₀ image of the
        # class-c contribution must equal the derived class-c2 contribution
        Ct = randn(rng, CT, b.mb, b.mo)
        Dc = sum(b.Bb[j] * Ct * b.Bo[j]' for j in 1:b.dim)
        Dc2 = sum(Bb2[j] * conj.(Ct) * Bo2[j]' for j in 1:b.dim)
        norm(da * conj.(Dc * Aa) - Dc2) < rtol * norm(Dc) || continue
        out[c] = SchurBlock{T}(b.dim, b.mb, b.mo, b.mf, b.Bb, b.Bo, Int8(1), c2)
        out[c2] = SchurBlock{T}(b2.dim, b2.mb, b2.mo, b2.mf, Bb2, Bo2, Int8(2), c)
    end
    return out
end

"""
    $(SIGNATURES)

Build the Schur-adapted bases for every IBZ kpoint. `frozen_ibz` is the
frozen-band mask at the IBZ kpoints. Errors when the frozen window cuts a
degenerate multiplet (the frozen subspace must be little-group invariant) or
when the feasibility counts `m_f ≤ m_o ≤ m_b` fail for some irrep.
"""
function schur_basis(
        sc::SymmetryConstraint{T}, frozen_ibz::AbstractMatrix{Bool}
    ) where {T}
    CT = Complex{T}
    rng = MersenneTwister(20260819)
    arng = MersenneTwister(20260820)   # separate stream: keep `rng` draws stable
    nb, nw = sc.nbands, sc.nwann
    blocks = Vector{Vector{SchurBlock{T}}}(undef, sc.nk_ibz)
    aop = Vector{Union{Nothing, Tuple{Matrix{CT}, Matrix{CT}}}}(undef, sc.nk_ibz)

    for iki in 1:sc.nk_ibz
        entries = sc.proj[iki]
        uidx = findall(e -> !e[3], entries)
        aidx = findall(e -> e[3], entries)
        aop[iki] = isempty(aidx) ? nothing :
            (Matrix(entries[aidx[1]][1]), Matrix(entries[aidx[1]][2]))

        # orbital representation W(ĥ) = A(ĥ)†
        Wo = [Matrix(entries[i][2])' for i in uidx]
        db = [Matrix(entries[i][1]) for i in uidx]

        ok = sc.band_ok[iki]
        f = frozen_ibz[:, iki] .& ok
        r = ok .& .!f
        # the frozen subspace must be invariant (energy blocks)
        for d in db
            norm(d[f, r]) < 1.0e-5 * max(1, norm(d)) ||
                error("frozen window cuts a degenerate multiplet at IBZ kpoint $iki")
        end

        # copies: orbital side (defines the references), then band side
        ocopies = _irrep_copies([Matrix{CT}(w) for w in Wo], rng)
        classes = Vector{Tuple{Vector{Matrix{CT}}, Vector{Matrix{CT}}, Vector{Matrix{CT}}, Int}}()
        # per class: (ρref, orbital copies Qo, band copies Qb, n frozen band copies)
        for (Q, ρQ) in ocopies
            matched = false
            for cl in classes
                u = _align_to(ρQ, cl[1], randn(rng, CT, size(Q, 2), size(Q, 2)))
                if u !== nothing
                    push!(cl[2], Q * u)
                    matched = true
                    break
                end
            end
            matched || push!(classes, (ρQ, [Q], Matrix{CT}[], 0))
        end
        # band side: frozen subspace first, then the rest
        for (sub, is_froz) in ((f, true), (r, false))
            count(sub) == 0 && continue
            idx = findall(sub)
            ρsub = [d[idx, idx] for d in db]
            for (Qs, ρQ) in _irrep_copies(ρsub, rng)
                Q = zeros(CT, nb, size(Qs, 2))
                Q[idx, :] .= Qs
                for (ci, cl) in enumerate(classes)
                    u = _align_to(ρQ, cl[1], randn(rng, CT, size(Qs, 2), size(Qs, 2)))
                    if u !== nothing
                        if is_froz
                            insert!(cl[3], cl[4] + 1, Q * u)
                            classes[ci] = (cl[1], cl[2], cl[3], cl[4] + 1)
                        else
                            push!(cl[3], Q * u)
                        end
                        break
                    end
                end
                # unmatched band copies carry irreps absent from the orbitals:
                # a covariant gauge has no weight there; drop them.
            end
        end

        blks = SchurBlock{T}[]
        for (ρref, Qos, Qbs, nf) in classes
            dλ = size(ρref[1], 1)
            mo, mb, mf = length(Qos), length(Qbs), nf
            mf <= mo <= mb || error(
                "infeasible symmetry constraint at IBZ kpoint $iki: " *
                    "irrep with (m_f, m_o, m_b) = ($mf, $mo, $mb) violates m_f ≤ m_o ≤ m_b"
            )
            Bo = [hcat((Q[:, j] for Q in Qos)...) for j in 1:dλ]
            Bb = [hcat((Q[:, j] for Q in Qbs)...) for j in 1:dλ]
            push!(blks, SchurBlock{T}(dλ, mb, mo, mf, Bb, Bo))
        end
        if aop[iki] !== nothing
            blks = _classify_aop(blks, aop[iki][1], aop[iki][2], arng)
            # drop the soft coset average where every class is exact
            all(b.akind != Int8(0) for b in blks) && (aop[iki] = nothing)
        end
        blocks[iki] = blks
    end

    nx = sum(sum(_block_nparams(blk)) for blks in blocks for blk in blks)
    return SchurBasis{T}(blocks, aop, nx)
end

# REAL parameter counts (nX, nY) of a block: complex entries as re/im pairs
# for unconstrained/pairing-source blocks, one real per entry for real-type
# blocks, 4 reals per quaternion (= one real per complex entry) for
# quaternionic blocks; derived pairing partners carry none
function _block_nparams(blk::SchurBlock)
    blk.akind == Int8(2) && return (0, 0)
    nX = blk.mo^2
    nY = (blk.mb - blk.mf) * (blk.mo - blk.mf)
    fac = blk.akind == Int8(3) || blk.akind == Int8(4) ? 1 : 2
    return (fac * nX, fac * nY)
end

# iterate the flat REAL parameter vector over the parameter-bearing blocks:
# yields (iki, blk, Xrange, Yrange)
function _schur_ranges(sb::SchurBasis)
    out = Tuple{Int, SchurBlock, UnitRange{Int}, UnitRange{Int}}[]
    off = 0
    for (iki, blks) in enumerate(sb.blocks), blk in blks
        blk.akind == Int8(2) && continue
        nX, nY = _block_nparams(blk)
        push!(out, (iki, blk, (off + 1):(off + nX), (off + nX + 1):(off + nX + nY)))
        off += nX + nY
    end
    return out
end

_blockC(blk, X, Y) = vcat(X[1:blk.mf, :], Y * X[(blk.mf + 1):end, :])

# materialize / store the (X, Y) block matrices from / into the flat REAL
# parameter vector, per akind storage layout. `grad = true` applies the
# gradient convention: real-type gradients are Re-projected, quaternionic
# component gradients carry the factor 2 from the two copies of each
# component in the structured matrix.
function _schur_getM(x::AbstractVector{T}, blk, r, m1::Int, m2::Int) where {T}
    if blk.akind == Int8(3)
        return Matrix(reshape(view(x, r), m1, m2))
    elseif blk.akind == Int8(4)
        v = reinterpret(Complex{T}, view(x, r))
        h = (m1 ÷ 2) * (m2 ÷ 2)
        return _quat_assemble(
            reshape(v[1:h], m1 ÷ 2, m2 ÷ 2), reshape(v[(h + 1):end], m1 ÷ 2, m2 ÷ 2)
        )
    else
        return Matrix(reshape(reinterpret(Complex{T}, view(x, r)), m1, m2))
    end
end

function _schur_setM!(x::AbstractVector{T}, blk, r, M; grad::Bool = false) where {T}
    if blk.akind == Int8(3)
        reshape(view(x, r), size(M)) .= real.(M)
    elseif blk.akind == Int8(4)
        A, B = _quat_extract(M)
        grad && (A .*= 2; B .*= 2)
        v = reinterpret(Complex{T}, view(x, r))
        h = length(A)
        v[1:h] .= vec(A)
        v[(h + 1):end] .= vec(B)
    else
        reshape(reinterpret(Complex{T}, view(x, r)), size(M)) .= M
    end
    return x
end

_schur_getX(x, blk, rX) = _schur_getM(x, blk, rX, blk.mo, blk.mo)
_schur_getY(x, blk, rY) = _schur_getM(x, blk, rY, blk.mb - blk.mf, blk.mo - blk.mf)

# anti-unitary coset average and its real-inner-product adjoint
_aavg(U, ::Nothing) = U
_aavg(U, da_Aa) = (U .+ da_Aa[1] * conj.(U * da_Aa[2])) ./ 2
_aavg_adj(G, ::Nothing) = G
_aavg_adj(G, da_Aa) = (G .+ transpose(da_Aa[1]) * conj.(G) * da_Aa[2]') ./ 2

"""
    $(SIGNATURES)

Decode the (real) Schur parameters `x` into the covariant IBZ gauge `U_ibz`.
"""
function schur_decode!(
        U_ibz::AbstractArray{<:Complex, 3}, x::AbstractVector{<:Real}, sb::SchurBasis
    )
    fill!(U_ibz, 0)
    for (iki, blk, rX, rY) in _schur_ranges(sb)
        X = _schur_getX(x, blk, rX)
        Y = _schur_getY(x, blk, rY)
        C = _blockC(blk, X, Y)
        for j in 1:blk.dim
            mul!(view(U_ibz, :, :, iki), blk.Bb[j] * C, blk.Bo[j]', true, true)
        end
        if blk.akind == Int8(1)
            pb = sb.blocks[iki][blk.partner]
            Cc = conj.(C)
            for j in 1:pb.dim
                mul!(view(U_ibz, :, :, iki), pb.Bb[j] * Cc, pb.Bo[j]', true, true)
            end
        end
    end
    for iki in axes(U_ibz, 3)
        view(U_ibz, :, :, iki) .= _aavg(U_ibz[:, :, iki], sb.aop[iki])
    end
    return U_ibz
end

"""
    $(SIGNATURES)

Chain the canonical IBZ gradient `G_ibz = dΩ/dU*` into the Schur parameter
gradient `g` (layout of `x`).
"""
function schur_encode_gradient!(
        g::AbstractVector{<:Real}, G_ibz::AbstractArray{<:Complex, 3},
        x::AbstractVector{<:Real}, sb::SchurBasis,
    )
    Ga = [_aavg_adj(G_ibz[:, :, iki], sb.aop[iki]) for iki in axes(G_ibz, 3)]
    for (iki, blk, rX, rY) in _schur_ranges(sb)
        X = _schur_getX(x, blk, rX)
        Y = _schur_getY(x, blk, rY)
        GC = zeros(eltype(G_ibz), blk.mb, blk.mo)
        for j in 1:blk.dim
            mul!(GC, blk.Bb[j]', Ga[iki] * blk.Bo[j], true, true)
        end
        if blk.akind == Int8(1)
            pb = sb.blocks[iki][blk.partner]
            GCp = zeros(eltype(G_ibz), blk.mb, blk.mo)
            for j in 1:pb.dim
                mul!(GCp, pb.Bb[j]', Ga[iki] * pb.Bo[j], true, true)
            end
            GC .+= conj.(GCp)
        end
        # C = J(Y) X: GX = J† GC ; GY = (GC X†) bottom-right block
        GX = vcat(GC[1:blk.mf, :], Y' * GC[(blk.mf + 1):end, :])
        GCXt = GC * X'
        GY = GCXt[(blk.mf + 1):end, (blk.mf + 1):end]
        _schur_setM!(g, blk, rX, GX; grad = true)
        _schur_setM!(g, blk, rY, GY; grad = true)
    end
    return g
end

"""
    $(SIGNATURES)

Initial Schur parameters from an IBZ gauge (its covariant block content).
"""
function schur_initial_x(
        U_ibz::AbstractArray{CT, 3}, sb::SchurBasis{T}
    ) where {CT <: Complex, T}
    x = zeros(T, sb.nx)
    for (iki, blk, rX, rY) in _schur_ranges(sb)
        C = zeros(CT, blk.mb, blk.mo)
        for j in 1:blk.dim
            mul!(C, blk.Bb[j]', view(U_ibz, :, :, iki) * blk.Bo[j], true, true)
        end
        C ./= blk.dim
        # direct (empty-shape-safe) construction of the per-block (X, Y):
        # Y spans the dominant non-frozen row space of C, X aligns J(Y)† C;
        # for constrained blocks C is first projected onto the structure and
        # Y is built structure-preserving (real / Kramers-paired columns)
        mf, mo, mb = blk.mf, blk.mo, blk.mb
        if blk.akind == Int8(3)
            Cs = real.(C)
            Y = zeros(T, mb - mf, mo - mf)
            if mo > mf && mb > mf
                Cr = Cs[(mf + 1):end, :]
                E, V = eigen(Symmetric(Cr * transpose(Cr)))
                Y .= V[:, (end - (mo - mf) + 1):end]
            end
        elseif blk.akind == Int8(4)
            Jb = cat(_jmat(CT, mf), _jmat(CT, mb - mf); dims = (1, 2))
            Cs = (C .+ Jb * conj.(C) * transpose(_jmat(CT, mo))) ./ 2
            Y = zeros(CT, mb - mf, mo - mf)
            if mo > mf && mb > mf
                Y .= _quat_range(Cs[(mf + 1):end, :], mo - mf)
            end
        else
            Cs = C
            Y = zeros(CT, mb - mf, mo - mf)
            if mo > mf && mb > mf
                Cr = Cs[(mf + 1):end, :]
                E, V = eigen(Hermitian(Cr * Cr'))
                Y .= V[:, (end - (mo - mf) + 1):end]
            end
        end
        JtC = vcat(Cs[1:mf, :], Y' * Cs[(mf + 1):end, :])
        X = orthonorm_lowdin(JtC)
        _schur_setM!(x, blk, rX, X)
        _schur_setM!(x, blk, rY, Y)
    end
    return x
end

"""Product-of-Stiefel manifold over the per-(kpoint, irrep) Schur blocks."""
struct SchurManifold <: Optim.Manifold
    ranges::Vector{Tuple{Int, SchurBlock, UnitRange{Int}, UnitRange{Int}}}
end
SchurManifold(sb::SchurBasis) = SchurManifold(_schur_ranges(sb))

function Optim.retract!(M::SchurManifold, x)
    for (_, blk, rX, rY) in M.ranges
        _schur_setM!(x, blk, rX, orthonorm_lowdin(_schur_getX(x, blk, rX)))
        if !isempty(rY)
            _schur_setM!(x, blk, rY, orthonorm_lowdin(_schur_getY(x, blk, rY)))
        end
    end
    return x
end

function Optim.project_tangent!(M::SchurManifold, g, x)
    for (_, blk, rX, rY) in M.ranges
        X = _schur_getX(x, blk, rX)
        GX = _schur_getX(g, blk, rX)
        GX .-= X * ((X' * GX .+ GX' * X) ./ 2)
        _schur_setM!(g, blk, rX, GX)
        if !isempty(rY)
            Y = _schur_getY(x, blk, rY)
            GY = _schur_getY(g, blk, rY)
            GY .-= Y * ((Y' * GY .+ GY' * Y) ./ 2)
            _schur_setM!(g, blk, rY, GY)
        end
    end
    return g
end
