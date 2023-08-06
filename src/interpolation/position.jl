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
    $(SIGNATURES)

Generate tight-binding position operator from a Wannierization [`Model`](@ref).

# Keyword Arguments
- `imlog_diag` and `force_hermiticity`: See [`compute_berry_connection_kspace`](@ref)
- others see the keyword args of [`generate_Rspace`](@ref)
"""
function TBPosition(
    model::Model,
    gauges::AbstractVector=model.gauges;
    imlog_diag::Bool=true,
    force_hermiticity::Bool=default_w90_berry_position_force_hermiticity(),
    kwargs...,
)
    Rspace = generate_Rspace(model; kwargs...)
    # Wannier-gauge position operator in kspace, WYSV Eq. 44
    Aᵂ = compute_berry_connection_kspace(model, gauges; imlog_diag, force_hermiticity)
    # Wannier-gauge position operator in Rspace, WYSV Eq. 43
    A_R = fourier(model.kpoints, Aᵂ, Rspace)
    bare_Rspace, bare_A_R = simplify(Rspace, A_R)
    return TBPosition(bare_Rspace, bare_A_R)
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
    eigenvalues, U = eigen(H_k)

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

    return eigenvalues, U, dH, D
end

"""Compute Wannier-gauge ``Rα * < m0 | H | nR >``, i.e., right-hand side of YWVS Eq. 38."""
@inline function compute_hamiltonian_gradient_Rspace(hamiltonian::TBOperator)
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
    return RH_R
end

@inline function compute_D_matrix(
    hamiltonian::TBOperator, kpoints::AbstractVector; kwargs...
)
    # k-space Hamiltonian
    H_k = invfourier(hamiltonian, kpoints)
    # Rα * < m0 | H | nR >
    RH_R = compute_hamiltonian_gradient_Rspace(hamiltonian)
    RH_k = invfourier(hamiltonian.Rspace, RH_R, kpoints)
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
    hamiltonian::TBOperator, kpoints::AbstractVector, fermi_energy::Real; kwargs...
)
    H_k = invfourier(hamiltonian, kpoints)
    # Rα * < m0 | H | nR >
    RH_R = compute_hamiltonian_gradient_Rspace(hamiltonian)
    RH_k = invfourier(hamiltonian.Rspace, RH_R, kpoints)
    return compute_J_matrix(H_k, RH_k, kpoints, fermi_energy; kwargs...)
end
