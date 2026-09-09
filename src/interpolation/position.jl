export TBPosition, position_transl_inv_full

"""
Construct a tight-binding position operator in R-space.

!!! note

    The tight-binding operator defined on the Rspace domain.
    Here the inner type `MVec3` represents 3 Cartesian directions.
    It is also possible to use `Vec3`, however, this will forbid all the
    in-place functions.
"""
function TBPosition end

function TBPosition(Rspace::BareRspace, operator::AbstractVector)
    @assert !isempty(operator) "empty operator"
    @assert !isempty(operator[1]) "empty operator"
    @assert operator[1] isa AbstractMatrix "operator must be a matrix"
    v = operator[1][1, 1]
    @assert v isa AbstractVector && length(v) == 3 "each element must be 3-vector"
    T = real(eltype(v))
    M = Matrix{MVec3{Complex{T}}}
    return TBOperator{M}("Position", Rspace, operator)
end

"""
    $(SIGNATURES)

Generate tight-binding position operator from a Wannierization [`Model`](@ref).

# Keyword Arguments
- `imlog_diag` and `force_hermiticity`: See [`compute_berry_connection_kspace`](@ref)
- `transl_inv_full`: use the translationally-invariant position operator,
    wannier90's `transl_inv_full`. See [`position_transl_inv_full`](@ref).
- others see the keyword args of [`generate_Rspace`](@ref)

!!! note

    When wannier90 writes a `tb.dat` file, it does not force the hermiticity of
    position operator, so to reproduce the same position operator in `tb.dat`,
    one should set `force_hermiticity=false`. However, when wannier90 computes
    Berry curvature, it does force the hermiticity when constructing position
    operator from `mmn` file.
"""
function TBPosition(
        Rspace::Union{WignerSeitzRspace, MDRSRspace},
        model::Model,
        gauges::AbstractVector = model.gauges;
        imlog_diag::Union{Nothing, Bool} = nothing,
        force_hermiticity::Union{Nothing, Bool} = nothing,
        transl_inv_full::Bool = false,
        kwargs...,
    )
    if transl_inv_full
        bare_Rspace, bare_A_R = position_transl_inv_full(
            Rspace, model, gauges; imlog_diag, force_hermiticity
        )
        return TBPosition(bare_Rspace, bare_A_R)
    end
    imlog_diag = something(imlog_diag, true)
    force_hermiticity = something(
        force_hermiticity, default_w90_berry_position_force_hermiticity()
    )
    # Wannier-gauge position operator in kspace, WYSV Eq. 44
    Aᵂ = compute_berry_connection_kspace(model, gauges; imlog_diag, force_hermiticity)
    # Wannier-gauge position operator in Rspace, WYSV Eq. 43
    A_R = fourier(model.kpoints, Aᵂ, Rspace)
    bare_Rspace, bare_A_R = simplify(Rspace, A_R)
    return TBPosition(bare_Rspace, bare_A_R)
end

"""
    $(SIGNATURES)

Translationally-invariant Wannier-gauge position operator in Rspace, i.e.
wannier90's `transl_inv_full`.

Returns a `(bare_Rspace, bare_A_R)` tuple.

The naive discretization of ``\\langle 0i | r | Rj \\rangle`` changes by more
than the shift itself when all Wannier centres are displaced by a constant
vector. Measuring the position relative to the midpoint
``\\bar{r}_{ij;R} = (r_i + r_j - R) / 2`` removes that error, at the cost of a
phase depending on both ``\\mathbf{b}`` and ``\\mathbf{R}``:

```math
A_{ij}(\\mathbf{R}) = \\sum_{\\mathbf{b}}
    e^{-i \\mathbf{b} \\cdot \\mathbf{R} / 2}
    \\mathcal{F}\\left[ e^{i \\mathbf{b} \\cdot (r_i + r_j) / 2}
    \\, i w_b \\mathbf{b} \\, M^{\\mathbf{k}, \\mathbf{b}}_{ij} \\right](\\mathbf{R})
```

Because the ``\\mathbf{R}``-dependent phase sits outside the Fourier transform,
the ``\\mathbf{b}``-sum can no longer be collapsed at each kpoint before
transforming; it becomes the outermost loop.

# Keyword Arguments
- `imlog_diag`: must be `false`. Rewriting the band-diagonal elements with
    ``\\mathrm{Im} \\ln M_{nn}`` (MV1997 Eq. (31), wannier90's `transl_inv`)
    is incompatible with this method: under a rigid shift of the Wannier
    centres by ``\\bm{t}`` the overlaps pick up
    ``M \\to e^{-i \\bm{b} \\cdot \\bm{t}} M``, which the midpoint phase
    cancels exactly, but the logarithm turns that factor into an additive
    ``+ \\bm{b} \\cdot \\bm{t}`` that no phase can cancel. wannier90 rejects
    the same combination in `postw90_readwrite.F90:857`. The band-diagonal at
    ``\\bm{R} = 0`` is supplied by the Wannier centres instead.
- `force_hermiticity`: must be `false`. postw90 hermitizes only on the
    non-translationally-invariant path (`get_oper.F90:743-749`), so
    hermitizing here would not reproduce `AA_R`.

!!! note

    The ``\\mathbf{b}``-sum is taken outside the Fourier transform, so slot `ib`
    must denote the same ``\\mathbf{b}``-vector at every kpoint. The stencil is
    therefore reordered to the Γ-point ordering by [`reorder`](@ref) and the
    overlaps are permuted to match, which is wannier90's `kmesh_info%nnord`.

    The ``e^{-i \\mathbf{b} \\cdot \\mathbf{R} / 2}`` phase is applied on the
    *simplified* R-vectors, i.e. after Wigner-Seitz degeneracy division and MDRS
    expansion. Those are postw90's `crvec_pw90`, built from
    `ws_distance%irdist`, i.e. the MDRS-shifted ``\\mathbf{R} + \\mathbf{T}``
    rather than the bare `irvec` (`postw90_common.F90:1832-1836, 1914-1917`);
    `operator_wigner_setup` applies the degeneracies before the phase, just as
    `reducer` does here. The result therefore agrees with postw90's `AA_R`.

    Both phases are transcribed from `get_oper.F90`: ``e^{+i b \\cdot \\bar{r}}``
    at `:661-664`, and ``e^{i (-1/2) b \\cdot \\mathrm{crvec}}`` at `:700-704`
    with `bk(:,nn,1)`, the first-kpoint b-vector, which is what `bvectors[ib]`
    is here. The translational-invariance test does not constrain the second
    phase — it contains neither ``\\bar{r}`` nor the shift — so its sign and
    its factor of one half rest on that transcription.
"""
function position_transl_inv_full(
        Rspace::Union{WignerSeitzRspace, MDRSRspace},
        model::Model,
        gauges::AbstractVector = model.gauges;
        imlog_diag::Union{Nothing, Bool} = nothing,
        force_hermiticity::Union{Nothing, Bool} = nothing,
    )
    # Both have exactly one legal value here, so imply them and reject only an
    # explicit request for the other one. wannier90 is in the same position:
    # `transl_inv` defaults to false, so `transl_inv_full = .true.` alone is a
    # valid input there.
    something(imlog_diag, false) && error(
        "transl_inv_full needs imlog_diag=false: the Im-log band-diagonal " *
        "breaks translational invariance, since it turns the overlap phase " *
        "into an additive shift that the midpoint phase cannot cancel"
    )
    something(force_hermiticity, false) && error(
        "transl_inv_full needs force_hermiticity=false: postw90 hermitizes " *
        "only on the non-translationally-invariant path " *
        "(get_oper.F90:743-749), so hermitizing here would not reproduce AA_R"
    )
    # slot `ib` must mean the same b-vector at every kpoint, since the b-sum is
    # moved outside the Fourier transform
    kstencil = reorder(model.kstencil)
    overlaps = reorder(
        model.overlaps, model.kstencil.kpb_k, model.kstencil.kpb_G, kstencil
    )

    kpoints = kstencil.kpoints
    nkpts = length(kpoints)
    nwann = size(gauges[1], 2)
    wb = kstencil.bweights
    C = complex(eltype(gauges[1]))

    # WF centres, MV1997 Eq. (31); wannier90's `wannier_centres_from_AA_R`
    centers = omega(kstencil, overlaps, gauges).r
    # midpoint r̄_ij = (r_i + r_j) / 2
    r̄ = [(centers[i] + centers[j]) / 2 for i in 1:nwann, j in 1:nwann]

    # the simplified R-vectors depend only on `Rspace`, so build the reducer
    # once and reuse it for every b-vector
    reducer = RvectorReducer(Rspace)
    bare_Rspace = BareRspace(Rspace.lattice, reducer.Rvectors)
    Rcarts = map(R -> bare_Rspace.lattice * R, bare_Rspace.Rvectors)

    bare_A_R = [zeros(Vec3{C}, nwann, nwann) for _ in 1:n_Rvectors(bare_Rspace)]

    for (ib, b) in enumerate(kstencil.bvectors)
        # `reorder` makes slot `ib` denote `b` at every kpoint
        # e^{i b ⋅ r̄_ij}, the phase making the b-sum invariant under a shift
        phase_r = map(v -> cis(b ⋅ v), r̄)

        Aᵂ = map(1:nkpts) do ik
            Uₖ = gauges[ik]
            Uₖ₂ = gauges[kstencil.kpb_k[ik][ib]]
            Mᵂ = Uₖ' * overlaps[ik][ib] * Uₖ₂
            # no `- I` here, unlike `compute_berry_connection_kspace`: wannier90
            # uses the overlap directly (get_oper.F90:612-613). In the collapsed
            # sum `∑_b w_b b = 0` kills the identity term, but here each b is
            # phased separately, so it would survive as a spurious R = 0 constant
            return wb[ib] .* Ref(b) .* (im .* Mᵂ .* phase_r)
        end

        A_R = reducer(fourier(kpoints, Aᵂ, Rspace))
        for (iR, Rcart) in enumerate(Rcarts)
            # e^{-i b ⋅ R / 2}, completing the midpoint (r_i + r_j - R) / 2
            bare_A_R[iR] .+= A_R[iR] .* cis(-(b ⋅ Rcart) / 2)
        end
    end

    # pin the R = 0 diagonal to the WF centres, so the position operator can
    # never disagree with the spread; wannier90 does the same
    iR0 = findfirst(iszero, bare_Rspace.Rvectors)
    isnothing(iR0) && error("R = 0 missing from the Rspace, cannot pin the centres")
    for i in 1:nwann
        bare_A_R[iR0][i, i] = centers[i]
    end

    return bare_Rspace, bare_A_R
end

function TBPosition(model::Model, gauges::AbstractVector = model.gauges; kwargs...)
    Rspace = generate_Rspace(model)
    return TBPosition(Rspace, model, gauges; kwargs...)
end

"""
    $(TYPEDEF)

A struct for interpolating tight-binding position operator on given kpoints.

# Fields
$(FIELDS)
"""
struct PositionInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian.
    Since we interpolate on kpoints in Bloch gauge, we need to store the Hamiltonain.
    """
    hamiltonian::TBOperator

    """R-space Hamiltonian gradient, Rα * < m0 | H | nR >, i.e., RHS of YWVS Eq. 38.
    Can be computed from hamiltonian operator by [`TBHamiltonianGradient`](@ref)."""
    hamiltonian_gradient::TBOperator

    """R-space position operator."""
    position::TBOperator
end

"""Interpolate the Hamiltonian operator and transform it to Bloch gauge."""
function (interp::PositionInterpolator)(
        kpoints::AbstractVector{<:AbstractVector}; kwargs...
    )
    _, gauges, _, D_matrices = compute_D_matrix(
        interp.hamiltonian, interp.hamiltonian_gradient, kpoints; kwargs...
    )

    # gauge-covariant part of k-space position operator
    Aᵂ_k = invfourier(interp.position, kpoints)
    # build the gauge-covariant position operator
    A_k = map(zip(Aᵂ_k, gauges, D_matrices)) do (Aᵂ, U, D)
        U' * Aᵂ * U + im * D
    end
    return A_k
end

"""
    $(SIGNATURES)

Compute the matrix D in YWVS Eq. 25 (or Eq. 32 if `degen_pert = true`).

# Arguments
- `kpoints`: fractional kpoints coordinates to be interpolated on

# Keyword arguments
- `degen_pert`: use perturbation treatment for degenerate eigenvalues
- `degen_tol`: degeneracy threshold in eV

# Return
- `eigenvalues`: energy eigenvalues
- `U`: the unitary transformation matrix
- `dH`: the covariant part of derivative of Hamiltonian in Bloch gauge,
    the ``\\bar{H}_{\\alpha}^{(H)}`` in YWVS Eq. 26
- `D`: the matrix ``D_{nm,\\alpha}^{(H)} = (U^\\dagger \\partial_{\\alpha}) U)_{nm},
    i.e., YWVS Eq. 25 or Eq. 32

!!! warning

    If `degen_pert = true`, the degenerate subspace is rotated such that
    ``\\bar{H}_{\\alpha}^{(H)}`` is diagonal, note only the ``\\alpha=x``
    direction is treated, since in general it is not possible to diagonalize
    simultaneously all the three directions.
"""
function compute_D_matrix(
        H_k::AbstractVector,
        RH_k::AbstractVector,
        kpoints::AbstractVector;
        degen_pert::Bool = default_w90_berry_use_degen_pert(),
        degen_tol::Real = default_w90_berry_degen_tol(),
    )
    nkpts = length(kpoints)
    @assert nkpts == length(H_k) == length(RH_k) > 0 "kpoints mismatched"
    nwann = size(H_k[1], 1)
    T = eltype(H_k[1])

    # first, need Hamiltonian eigenvalues and eigenvectors
    eigenvalues, U = eigen(H_k)

    # the covariant part of Hamiltonian gauge dH, i.e., dHᴴ = U† dHᵂ U
    # also the ``\bar{H}_{\alpha}^{(H)}`` in YWVS Eq. 26
    # inner MVec3 for the three Cartesian directions
    dH = [zeros(MVec3{T}, nwann, nwann) for _ in 1:nkpts]
    # the D matrix = U† ∂U in Hamiltonian gauge, i.e. YWVS Eq. 25 or Eq. 32
    D = [zeros(MVec3{T}, nwann, nwann) for _ in 1:nkpts]

    for ik in 1:nkpts
        # derivative of Hamiltonian dH = [dH/dkx, dH/dky, dH/dkz]
        # in Wannier gauge, at kpoint k
        dHᵂₖ = RH_k[ik]
        # to Bloch gauge, dHₖ = U† dHᵂₖ U
        Uₖ = U[ik]
        dH[ik] .= Uₖ' * dHᵂₖ * Uₖ

        # the D matrix
        Δε = eigenvalues[ik] .- eigenvalues[ik]'
        # assign a nonzero number to the diagonal elements for inversion
        Δε[diagind(Δε)] .= 1
        Dₖ = D[ik]
        Dₖ .= dH[ik] ./ (-Δε)
        Dₖ[diagind(Dₖ)] .= Ref([0, 0, 0])

        # TODO: maybe it is helpful to run at least once the perturbation treatment
        # for one Cartesian direction, to avoid vanishing denominator in D matrix
        degen_pert || continue

        # now considering possible degeneracies
        # eigenvalues[ik][mask] are eigenvalues to be checked
        mask = trues(nwann)
        while any(mask)
            e = eigenvalues[ik][mask][1]
            # indices of degenerate eigenvalues
            idx = abs.(eigenvalues[ik] .- e) .< degen_tol
            if count(idx) > 1
                # I can only run once the diagonalization for only one Cartesian
                # direction, and update the U matrix. The following directions
                # will use the updated U matrix, and I only set the D matrix to
                # zero for the degenerate subspace.
                α = 1
                # diagonalize the submatrix
                h = map(x -> x[α], dH[ik][idx, idx])
                v, u = eigen(h)
                # update U such that in Hamiltonian gauge both H and dH
                # are diagonal in the degenerate subspace
                U[ik][idx, idx] *= u
                for (i0, i) in enumerate((1:nwann)[idx])
                    for j in (1:nwann)[idx]
                        if i == j
                            dH[ik][i, j][α] = v[i0]
                        else
                            dH[ik][i, j][α] = 0
                        end
                        # the D matrix
                        D[ik][i, j][α] = 0
                    end
                end
            end
            mask[idx] .= false
        end
    end

    return eigenvalues, U, dH, D
end

"""Compute Wannier-gauge ``Rα * < m0 | H | nR >``, i.e., right-hand side of YWVS Eq. 38."""
@inline function TBHamiltonianGradient(hamiltonian::TBOperator)
    # R-space Hamiltonain
    H_R = hamiltonian
    lattice = real_lattice(H_R)
    RH_R = map(zip(H_R.Rspace, H_R)) do (R, H)
        # to Cartesian in angstrom, result indexed by RH_R[iR][m, n][α], where
        # - iR is the index of Rvectors
        # - m, n are the indices of Wannier functions
        # - α ∈ {1, 2, 3} is the Cartesian direction for x, y, z
        Ref(im * (lattice * R)) .* H
    end
    return TBOperator("HamiltonianGradient", H_R.Rspace, RH_R)
end

@inline function compute_D_matrix(
        hamiltonian::TBOperator,
        hamiltonian_gradient::TBOperator,
        kpoints::AbstractVector;
        kwargs...,
    )
    # k-space Hamiltonian
    H_k = invfourier(hamiltonian, kpoints)
    # Rα * < m0 | H | nR >
    RH_k = invfourier(hamiltonian_gradient, kpoints)
    return compute_D_matrix(H_k, RH_k, kpoints; kwargs...)
end

"""
Compute `J`, `J⁻` and `J⁺` matrices, i.e., LVTS12 Eq. 76 and 77.

The `Jᴴ` matrix (LVTS12 Eq. 75) `= im * Dᴴ` matrix. Here we compute the
`Jᴴ⁺` and `Jᴴ⁻` matrices, which take occupations into account so they are
numerically more stable than the `Jᴴ` matrix. Then the `J`, `J⁻` and `J⁺`
are the Wannier-gauge matrices rotated from the `Jᴴ`, `Jᴴ⁻` and `Jᴴ⁺` matrices.
"""
function compute_J_matrix end

"""
    $(SIGNATURES)

Compute `J` matrices from `D` matrices.

# Arguments
- `eigenvalues`, `U`, `Dᴴ`: return values of [`compute_D_matrix`](@ref)
"""
function compute_J_matrix(
        eigenvalues::AbstractVector{<:AbstractVector},
        U::AbstractVector{<:AbstractMatrix},
        Dᴴ::AbstractVector{<:AbstractMatrix},
        fermi_energy::Real,
    )
    # occupations in Hamiltonian gauge
    fᴴ = [Diagonal(Int.(εₖ .<= fermi_energy)) for εₖ in eigenvalues]
    gᴴ = Ref(I) .- fᴴ

    Jᴴ = im * Dᴴ
    Jᴴ⁻ = fᴴ .* Jᴴ .* gᴴ
    Jᴴ⁺ = gᴴ .* Jᴴ .* fᴴ

    Ut = adjoint.(U)
    J = U .* Jᴴ .* Ut
    J⁻ = U .* Jᴴ⁻ .* Ut
    J⁺ = U .* Jᴴ⁺ .* Ut
    return J, J⁻, J⁺
end

"""
    $(SIGNATURES)

# Keyword Arguments
See [`compute_D_matrix`](@ref).
"""
@inline function compute_J_matrix(
        H_k::AbstractVector{<:AbstractMatrix},
        RH_k::AbstractVector{<:AbstractMatrix},
        kpoints::AbstractVector{<:AbstractVector},
        fermi_energy::Real;
        kwargs...,
    )
    eigenvalues, U, _, Dᴴ = compute_D_matrix(H_k, RH_k, kpoints; kwargs...)
    return compute_J_matrix(eigenvalues, U, Dᴴ, fermi_energy)
end

@inline function compute_J_matrix(
        hamiltonian::TBOperator,
        hamiltonian_gradient::TBOperator,
        kpoints::AbstractVector,
        fermi_energy::Real;
        kwargs...,
    )
    H_k = invfourier(hamiltonian, kpoints)
    # Rα * < m0 | H | nR >
    RH_k = invfourier(hamiltonian_gradient, kpoints)
    return compute_J_matrix(H_k, RH_k, kpoints, fermi_energy; kwargs...)
end
