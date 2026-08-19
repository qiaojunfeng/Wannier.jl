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
        centers::AbstractVector,
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
            end
        end
    end
    any(iszero, ikb) && error("Some (ibi, iki) pass-0 entries were never filled")

    return SymmetryConstraint{T}(
        nk_fbz, nk_ibz, nbvecs, nwann, nbands,
        fbz2ibz, ibz2fbz, stars,
        proj, band_ok, proj_niter,
        Lmat, trev_f,
        ikb, Aib, trev_ib,
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
    # decode (X,Y) -> covariant U at IBZ -> full mesh
    XY_to_X_Y!(ws.X_ibz, ws.Y_ibz, xy)
    X_Y_to_U!(ws.U_ibz, ws.X_ibz, ws.Y_ibz)
    project_covariant!(ws.U_ibz, sc)
    expand_gauges!(ws.full.U, ws.U_ibz, sc)

    kstencil, overlaps = model.kstencil, model.overlaps
    compute_MU_UtMU!(ws.full, kstencil, overlaps, ws.full.U)

    if G !== nothing
        omega_grad!(ws.full.GU, ws.full, kstencil, overlaps)
        pullback_gauges!(ws.G_ibz, ws.full.GU, sc)
        project_covariant!(ws.G_ibz, sc)
        encode_gradient_xy!(G, ws.G_ibz, ws.X_ibz, ws.Y_ibz, ws.frozen_ibz)
    end

    F === nothing && return nothing
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
        zeros(CT, nw, nbv, nkf), zeros(Vec3{T}, nw),
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
    nw, nki, nkf, nbv = sc.nwann, sc.nk_ibz, sc.nk_fbz, sc.nbvecs
    wb = sc.bweights

    # decode -> covariant U at IBZ
    XY_to_X_Y!(ws.X_ibz, ws.Y_ibz, xy)
    X_Y_to_U!(ws.U_ibz, ws.X_ibz, ws.Y_ibz)
    project_covariant!(ws.U_ibz, sc)

    # Pass 0 (heavy, IBZ only): U(ki+bi), MU_i, M̃_i
    for iki in 1:nki, ibi in 1:nbv
        ikb = sc.ikb[ibi, iki]
        Uk = view(ws.Ukb, :, :, ibi, iki)
        mul!(Uk, view(ws.U_ibz, :, :, ikb), sc.Aib[ibi, iki])
        sc.trev_ib[ibi, iki] && (Uk .= conj.(Uk))
        mul!(view(ws.MU, :, :, ibi, iki), view(M_ibz, :, :, ibi, iki), Uk)
        mul!(
            view(ws.Mt, :, :, ibi, iki), view(ws.U_ibz, :, :, iki)',
            view(ws.MU, :, :, ibi, iki),
        )
    end

    # Sweep 1 (light, full mesh): transported diagonals and centers
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
                ws.r[n] -= imaglog(t) * (wb[ibf] * sc.bvec_cart[ibf])
            end
        end
    end
    ws.r ./= nkf

    Ω = nothing
    if F !== nothing
        Ω = zero(T)
        @inbounds for ikf in 1:nkf, ibf in 1:nbv, n in 1:nw
            t = ws.tdiag[n, ibf, ikf]
            Ω += wb[ibf] * (1 - abs2(t) + imaglog(t)^2)
        end
        Ω = Ω / nkf - sum(r -> sum(abs2, r), ws.r)
    end

    if G !== nothing
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
                    q = imaglog(t) + sc.bvec_cart[ibf] ⋅ ws.r[n]
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
        # assembly: dΩ/dU*(ki) = Σ_bi MU_i 𝒦, then project and encode
        fill!(ws.G_ibz, 0)
        for iki in 1:nki, ibi in 1:nbv
            mul!(
                view(ws.G_ibz, :, :, iki), view(ws.MU, :, :, ibi, iki),
                view(ws.K, :, :, ibi, iki), true, true,
            )
        end
        project_covariant!(ws.G_ibz, sc)
        encode_gradient_xy!(G, ws.G_ibz, ws.X_ibz, ws.Y_ibz, ws.frozen_ibz)
    end

    return Ω
end

# -----------------------------------------------------------------------------
# Optimization driver
# -----------------------------------------------------------------------------

export localize_symmetric

"""
    $(SIGNATURES)

Minimize the MV spread over symmetry-covariant gauges parameterized at the
IBZ kpoints only (the SAWF constrained problem). Returns `(U_fbz, U_ibz)`:
the optimized covariant gauge expanded to the full mesh, and its IBZ
representative.

# Arguments
- `model`: full-mesh model in the *global* b ordering (see
  [`globalize_stencil`](@ref)), with overlaps unfolded from the IBZ. Its
  `gauges` provide the starting point (their IBZ slices are projected onto
  the covariant subspace).
- `M_ibz`: IBZ overlaps (`.immn`), needed for `level = 2`.
- `sc`: the [`SymmetryConstraint`](@ref).

# Keyword arguments
- `level`: `2` (default) evaluates value/gradient via the IBZ-only transport
  kernels ([`symmetric_fg2!`](@ref)); `1` expands to the full mesh each
  iteration ([`symmetric_fg1!`](@ref)). Identical results, different cost.
- remaining kwargs are forwarded to [`OptimLBFGS`](@ref).
"""
function localize_symmetric(
        model::Model,
        M_ibz::AbstractArray{<:Complex, 4},
        sc::SymmetryConstraint;
        level::Integer = 2,
        kwargs...,
    )
    solver = OptimLBFGS(; kwargs...)

    U0_ibz = extract_ibz_gauges(model.gauges, sc)
    project_covariant!(U0_ibz, sc)
    frozen_ibz = model.frozen_bands[:, sc.ibz2fbz]
    X0, Y0 = U_to_X_Y(U0_ibz, frozen_ibz)
    x0 = X_Y_to_XY(X0, Y0)

    if level == 1
        ws1 = SymmetricWorkspace(model, sc)
        fg! = (F, G, x) -> symmetric_fg1!(F, G, x, model, sc, ws1)
    elseif level == 2
        ws2 = SymmetricWorkspace2(model.eigenvalues, model.frozen_bands, sc)
        fg! = (F, G, x) -> symmetric_fg2!(F, G, x, M_ibz, sc, ws2)
    else
        error("level must be 1 or 2")
    end

    nw, nb = sc.nwann, sc.nbands
    per_k = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw)
    )
    man = Optim.PowerManifold(per_k, (nw^2 + nb * nw,), (sc.nk_ibz,))

    opt = _run_optim_fg!(fg!, x0, man, solver)
    xy = Optim.minimizer(opt)

    X, Y = XY_to_X_Y(xy, nb, nw)
    U_ibz = X_Y_to_U(X, Y)
    project_covariant!(U_ibz, sc)
    U_fbz = expand_gauges(U_ibz, sc)
    return U_fbz, U_ibz
end
