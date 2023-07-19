export PositionRspace, PositionKspace

"""
    $(TYPEDEF)

A struct representing tight-binding position operator in R-space.

# Fields
$(FIELDS)
"""
struct PositionRspace{M<:AbstractMatrix} <: AbstractOperatorRspace
    """the R-space domain (or called R-vectors) on which the operator is defined"""
    domain::BareRspaceDomain

    """The tight-binding operator defined on the domain.
    Here the type `M` is often `Matrix{MVec3{ComplexF64}}`, where the inner
    `MVec3` represents 3 Cartesian directions. One can also use `Vec3`,
    however, this will prohibit all the in-place functions."""
    operator::Vector{M}
end

"""
    $(TYPEDEF)

A struct representing tight-binding position operator on a kpoint grid.

# Fields
$(FIELDS)
"""
struct PositionKspace{K<:AbstractKpointContainer,M<:AbstractMatrix} <:
       AbstractOperatorKspace{K}
    """a [`KpointGrid`](@ref) or [`KpointList`](@ref) on which the operator is defined"""
    domain::K

    """the tight-binding operator defined on the domain"""
    operator::Vector{M}
end

function transform_gauge(
    position_k::PositionKspace, gauges::AbstractVector, D_matrices::AbstractVector
)
    A = map(zip(position_k.operator, gauges, D_matrices)) do (Aᵂ, U, D)
        U' * Aᵂ * U + im * D
    end
    return HamiltonianKspace(position_k.domain, A)
end

function transform_gauge(position_k::PositionKspace, hamiltonian_R::HamiltonianRspace)
    _, U, _, D = compute_D_matrix(hamiltonian_R, position_k.domain)
    return transform_gauge(position_k, U, D)
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
    hamiltonian_R::HamiltonianRspace,
    kpoints::Union{AbstractVector,AbstractKpointContainer};
    degen_pert::Bool=default_w90_berry_use_degen_pert(),
    degen_tol::Real=default_w90_berry_degen_tol(),
)
    H_k = invfourier(hamiltonian_R, kpoints)
    RH_R = zeros(hamiltonian_R)
    lattice = real_lattice(hamiltonian_R)
    map(1:nRvecs) do iR
        # to Cartesian in angstrom
        RH_R[iR] .= im * (lattice * hamiltonian_R.domain[iR]) * hamiltonian_R[iR]
    end
    RH_k = invfourier(RH_R, kpoints)
    return compute_D_matrix(H_k, RH_k, kpoints; degen_pert, degen_tol)
end

function compute_D_matrix(
    H_k::AbstractVector,
    RH_k::AbstractVector,
    kpoints;
    degen_pert::Bool=default_w90_berry_use_degen_pert(),
    degen_tol::Real=default_w90_berry_degen_tol(),
)
    nkpts = length(kpoints)
    @assert nkpts == length(H_k) == length(RH_k) > 0 "kpoints mismatched"
    nwann = size(H_k[1], 1)
    T = eltype(H_k[1])

    # first, need Hamiltonian eigenvalues and eigenvectors
    eigvalvecs = eigen.(H_k)
    eigvals = [x.values for x in eigvalvecs]
    U = [x.vectors for x in eigvalvecs]

    # the covariant part of Hamiltonian gauge dH, i.e., dHᴴ = U† dHᵂ U
    # also the ``\bar{H}_{\alpha}^{(H)}`` in YWVS Eq. 26
    # inner MVec3 for the three Cartesian directions
    dH = [zeros(MVec3{T}, nwann, nwann) for _ in 1:nkpts]
    # the D matrix = U† ∂U in Hamiltonian gauge, i.e. YWVS Eq. 25 or Eq. 32
    D = [zeros(MVec3{T}, nwann, nwann) for _ in 1:nkpts]

    for ik in 1:nkpts
        # derivative of Hamiltonain dH = [dH/dkx, dH/dky, dH/dkz]
        # in Wannier gauge, at kpoint k
        dHᵂₖ = @view RH_k[ik]
        # to Bloch gauge, dHₖ = U† dHᵂₖ U
        Uₖ = @view U[ik]
        dH[ik] .= Uₖ' * dHᵂₖ * Uₖ

        # the D matrix
        Δε = eigvals[ik] .- eigvals[ik]'
        # assign a nonzero number to the diagonal elements for inversion
        Δε[diagind(Δε)] .= 1
        Dₖ = @view D[ik]
        Dₖ .= dH[ik] ./ (-Δε)
        Dₖ[diagind(Dₖ)] .= 0

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
