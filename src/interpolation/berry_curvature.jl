export BerryCurvatureInterpolator

"""
    $(TYPEDEF)

A struct for interpolating Berry curvature on given kpoints.

# Fields
$(FIELDS)
"""
struct BerryCurvatureInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian."""
    hamiltonian::TBOperator

    """R-space position operator."""
    position::TBOperator
end

"""Interpolate Berry curvature and transform it to Bloch gauge, WYSV Eq. 27."""
function (interp::BerryCurvatureInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}; kwargs...
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    _, U, _, Dᴴ = compute_D_matrix(interp.hamiltonian, kpoints; kwargs...)

    # WYSV Eq. 27
    # gauge-covariant part of k-space position operator
    Aᵂ = invfourier(interp.position, kpoints)
    Āᴴ = map(zip(Aᵂ, U)) do (Aᵂₖ, Uₖ)
        Uₖ' * Aᵂₖ * Uₖ
    end

    # commutators [Dᴴα, Āᴴβ] - [Dᴴβ, Āᴴα]
    Dᴴ_Āᴴ = map(zip(Dᴴ, Āᴴ)) do (Dᴴₖ, Āᴴₖ)
        # each element: 3x1 vector *  1x3 vector -> 3x3 matrix,
        # the result can be indexed by [m, n][α, β], where
        # - m, n are band indices
        # - α, β are Cartesian indices ∈ {1, 2, 3}
        t_A = transpose.(Āᴴₖ)
        t_D = transpose.(Dᴴₖ)
        DA = Dᴴₖ * t_A - transpose.(Āᴴₖ * t_D)
        return DA - transpose.(DA)
    end

    # commutator [Dᴴα, Dᴴβ]
    Dᴴ_Dᴴ = map(Dᴴ) do Dᴴₖ
        # the result can be indexed by [m, n][α, β]
        DD = Dᴴₖ * transpose.(Dᴴₖ)
        return DD - transpose.(DD)
    end

    Rr_R = compute_berry_curvature_Rspace(interp.hamiltonian, interp.position)
    Ωᵂ = invfourier(interp.hamiltonian.Rspace, Rr_R, kpoints)
    # gauge-covariant part of Berry curvature, Ω̄ᴴ
    Ω̄ᴴ = map(zip(Ωᵂ, U)) do (Ωᵂₖ, Uₖ)
        Uₖ' * Ωᵂₖ * Uₖ
    end

    Ωᴴ = Ω̄ᴴ .- Dᴴ_Āᴴ .- im .* Dᴴ_Dᴴ
    # just need band-diagonal elements
    return map(Ωᴴ) do Ωᴴₖ
        diag(Ωᴴₖ)
    end
end

"""
Wannier-gauge Berry curvature Ωᵂ in Rspace.

Right-hand side of WYSV Eq. 40, `i * Rα * <0n|r̂β|Rm> - i * Rβ * <0n|r̂α|Rm>`
"""
function compute_berry_curvature_Rspace(hamiltonian::TBOperator, position::TBOperator)
    lattice = real_lattice(position)
    r_R = position.operator
    # Rr_R indexed by [iR][m, n][α, β]
    Rr_R = map(zip(hamiltonian.Rspace, r_R)) do (R, r)
        # to Cartesian in angstrom
        Rr = Ref(im * (lattice * R)) .* transpose.(r)
        return Rr - transpose.(Rr)
    end
    return Rr_R
end

"""Interpolate Berry curvature and transform it to Bloch gauge, sumed over
bands with given occupations, WYSV Eq. 32.

More accurate than directly summing over the band-resolved Berry curvature
times occupation, since the occupation factors are reorganized so that
pairs of fully occupied states give zero contribution.
"""
function (interp::BerryCurvatureInterpolator)(
    kpoints::AbstractVector{<:AbstractVector},
    occupations::AbstractVector{<:AbstractVector};
    kwargs...,
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    nkpts = length(kpoints)
    @assert length(occupations) == nkpts > 0 "length of occupations != length of kpoints"
    @assert length(occupations[1]) == n_wannier(interp.hamiltonian) "n_wannier of occupations != n_wannier of Hamiltonian"
    _, U, _, Dᴴ = compute_D_matrix(interp.hamiltonian, kpoints; kwargs...)

    # WYSV Eq. 27
    # gauge-covariant part of k-space position operator
    Aᵂ = invfourier(interp.position, kpoints)
    Āᴴ = map(zip(Aᵂ, U)) do (Aᵂₖ, Uₖ)
        Uₖ' * Aᵂₖ * Uₖ
    end

    Rr_R = compute_berry_curvature_Rspace(interp.hamiltonian, interp.position)
    Ωᵂ = invfourier(interp.hamiltonian.Rspace, Rr_R, kpoints)
    # gauge-covariant part of Berry curvature, Ω̄ᴴ
    Ω̄ᴴ = map(zip(Ωᵂ, U)) do (Ωᵂₖ, Uₖ)
        Uₖ' * Ωᵂₖ * Uₖ
    end

    return map(zip(Ω̄ᴴ, Āᴴ, Dᴴ, occupations)) do (Ω̄ᴴₖ, Āᴴₖ, Dᴴₖ, fₖ)
        # occupation -> matrix
        F = transpose(fₖ) .- fₖ
        DA = (F .* Dᴴₖ) * transpose.(Āᴴₖ)
        return sum(fₖ .* diag(Ω̄ᴴₖ)) +
               sum(diag(DA - transpose.(DA) + (im .* F .* Dᴴₖ) * transpose.(Dᴴₖ)))
    end
end

"""
    $(SIGNATURES)

Convert a Berry curvature vector to second-rank antisymmetric tensor.

WYSV Eq. 5, Ωγ = ϵαβγ * Ωαβ

# Arguments
- `Ω`: Berry curvature vector, `Ω = (Ωx, Ωy, Ωz)`
"""
function vector_to_tensor(Ω::AbstractArray{<:Number})
    Ωx, Ωy, Ωz = Ω
    return MMat3([
        0 Ωz -Ωy
        -Ωz 0 Ωx
        Ωy -Ωx 0
    ])
end
