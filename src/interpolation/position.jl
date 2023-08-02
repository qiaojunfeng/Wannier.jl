export TBPosition

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

    """R-space position operator."""
    position::TBOperator
end

"""Interpolate the Hamiltonian operator and transform it to Bloch gauge."""
function (interp::PositionInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}; kwargs...
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    _, gauges, _, D_matrices = compute_D_matrix(interp.hamiltonian, kpoints; kwargs...)

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
    degen_pert::Bool=default_w90_berry_use_degen_pert(),
    degen_tol::Real=default_w90_berry_degen_tol(),
)
    nkpts = length(kpoints)
    @assert nkpts == length(H_k) == length(RH_k) > 0 "kpoints mismatched"
    nwann = size(H_k[1], 1)
    T = eltype(H_k[1])

    # first, need Hamiltonian eigenvalues and eigenvectors
    eigvals, U = eigen(H_k)

    # the covariant part of Hamiltonian gauge dH, i.e., dHᴴ = U† dHᵂ U
    # also the ``\bar{H}_{\alpha}^{(H)}`` in YWVS Eq. 26
    # inner MVec3 for the three Cartesian directions
    dH = [zeros(MVec3{T}, nwann, nwann) for _ in 1:nkpts]
    # the D matrix = U† ∂U in Hamiltonian gauge, i.e. YWVS Eq. 25 or Eq. 32
    D = [zeros(MVec3{T}, nwann, nwann) for _ in 1:nkpts]

    for ik in 1:nkpts
        # derivative of Hamiltonain dH = [dH/dkx, dH/dky, dH/dkz]
        # in Wannier gauge, at kpoint k
        dHᵂₖ = RH_k[ik]
        # to Bloch gauge, dHₖ = U† dHᵂₖ U
        Uₖ = U[ik]
        dH[ik] .= Uₖ' * dHᵂₖ * Uₖ

        # the D matrix
        Δε = eigvals[ik] .- eigvals[ik]'
        # assign a nonzero number to the diagonal elements for inversion
        Δε[diagind(Δε)] .= 1
        Dₖ = D[ik]
        Dₖ .= dH[ik] ./ (-Δε)
        Dₖ[diagind(Dₖ)] .= Ref([0, 0, 0])

        # TODO: maybe it is helpful to run at least once the perturbation treatment
        # for one Cartesian direction, to avoid vanishing denominator in D matrix
        degen_pert || continue

        # now considering possible degeneracies
        # eigvals[ik][mask] are eigenvalues to be checked
        mask = trues(nwann)
        while any(mask)
            e = eigvals[ik][mask][1]
            # indices of degenerate eigenvalues
            idx = abs.(eigvals[ik] .- e) .< degen_tol
            if count(idx) > 1
                # I can only run once the diagonalization for only one Cartesian
                # direction, and update the U matrix. The following directions
                # will use the updated U matrix, and I only set the D matrix to
                # zero for the degenerate subspace.
                α = 1
                # diagonalize the submatrix
                h = map(x -> x[α], dH[ik][idx, idx])
                v, u = eigen(h)
                # update U such that in Hamiltonain gauge both H and dH
                # are diagonal in the degenerate subspace
                U[ik][idx, idx] *= u
                for i in idx
                    for j in idx
                        if i == j
                            dH[ik][i, j][α] .= v[i]
                        else
                            dH[ik][i, j][α] .= 0
                        end
                        # the D matrix
                        D[ik][i, j][α] .= 0
                    end
                end
            end
            mask[idx] .= false
        end
    end

    return eigvals, U, dH, D
end

function compute_D_matrix(hamiltonian::TBOperator, kpoints::AbstractVector; kwargs...)
    # R-space Hamiltonain
    H_R = hamiltonian
    # k-space Hamiltonian
    H_k = invfourier(H_R, kpoints)

    RH_R = zeros(H_R)
    lattice = real_lattice(H_R)
    RH_R = map(zip(H_R.Rspace, H_R)) do (R, H)
        # to Cartesian in angstrom, result indexed by RH_R[iR][m, n][α], where
        # - iR is the index of Rvectors
        # - m, n are the indices of Wannier functions
        # - α ∈ {1, 2, 3} is the Cartesian direction for x, y, z
        Ref(im * (lattice * R)) .* H
    end
    RH_k = invfourier(hamiltonian.Rspace, RH_R, kpoints)

    return compute_D_matrix(H_k, RH_k, kpoints; kwargs...)
end
