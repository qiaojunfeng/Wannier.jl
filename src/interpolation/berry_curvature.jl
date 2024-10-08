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

    """R-space Hamiltonian gradient, Rα * < m0 | H | nR >, i.e., RHS of YWVS Eq. 38.
    Can be computed from hamiltonian operator by [`TBHamiltonianGradient`](@ref)."""
    hamiltonian_gradient::TBOperator

    """R-space position operator."""
    position::TBOperator

    """Wannier-gauge Berry curvature, WYSV06 Eq. 40 or LVTS12 Eq. 31.
    Can be computed from position operator by [`TBBerryCurvature`](@ref)."""
    berry_curvature::TBOperator

    """Fermi energy in eV"""
    fermi_energy::Float64
end

function BerryCurvatureInterpolator(
    hamiltonian::TBOperator, position::TBOperator, fermi_energy::Real
)
    hamiltonian_gradient = TBHamiltonianGradient(hamiltonian)
    berry_curvature = TBBerryCurvature(position)
    return BerryCurvatureInterpolator(
        hamiltonian, hamiltonian_gradient, position, berry_curvature, fermi_energy
    )
end

"""
Wannier-gauge Berry curvature Ωᵂ in Rspace.

Right-hand side of WYSV Eq. 40, `i * Rα * <0n|r̂β|Rm> - i * Rβ * <0n|r̂α|Rm>`
"""
function TBBerryCurvature(position::TBOperator)
    lattice = real_lattice(position)
    r_R = position.operator
    # Rr_R indexed by [iR][m, n][α, β]
    Rr_R = map(zip(position.Rspace, r_R)) do (R, r)
        # to Cartesian in angstrom
        Rr = Ref(im * (lattice * R)) .* transpose.(r)
        return Rr - transpose.(Rr)
    end
    v = position.operator[1][1, 1]
    T = real(eltype(v))
    M = Matrix{MMat3{Complex{T}}}
    return TBOperator{M}("BerryCurvature", position.Rspace, Rr_R)
end

abstract type AbstractBerryCurvatureInterpolationAlgorithm end

"""
Interpolate Berry curvature using Eq. 32 in
X. Wang, J. R. Yates, I. Souza, and D. Vanderbilt, Phys. Rev. B 74, 195118 (2006).

This algorithm returns Berry curvature summed over occupied states.
It is more accurate than manually summing over the band-resolved Berry curvature
(see [`WYSV06BandResolved`](@ref)) times occupation, since the occupation factors
are reorganized so that pairs of fully occupied states give zero contribution.
"""
struct WYSV06 <: AbstractBerryCurvatureInterpolationAlgorithm end

"""
Interpolate Berry curvature using Eq. 27 in
X. Wang, J. R. Yates, I. Souza, and D. Vanderbilt, Phys. Rev. B 74, 195118 (2006).

This algorithm returns band-resolved Berry curvature, which might be useful for
plotting Berry-curvature-colored band structure.
See also [`WYSV06`](@ref).
"""
struct WYSV06BandResolved <: AbstractBerryCurvatureInterpolationAlgorithm end

"""
Interpolate Berry curvature using Eq. 51 in
M. Lopez, D. Vanderbilt, T. Thonhauser, and I. Souza, Phys. Rev. B 85, 014435 (2012).

Compared to [`WYSV06`](@ref), this algorithm computes Berry curvature of the
whole occuppied manifold in Wannier gauge. The algorithm is probably more
numerically stable than [`WYSV06`](@ref). However, it does not return
band-resolved Berry curvature (see [`WYSV06BandResolved`](@ref)).
"""
struct LVTS12 <: AbstractBerryCurvatureInterpolationAlgorithm end

"""Interpolate Berry curvature and transform it to Bloch gauge, WYSV Eq. 27."""
function (interp::BerryCurvatureInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}, ::WYSV06BandResolved; kwargs...
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    _, U, _, Dᴴ = compute_D_matrix(
        interp.hamiltonian, interp.hamiltonian_gradient, kpoints; kwargs...
    )

    # WYSV Eq. 27
    # gauge-covariant part of k-space position operator
    Aᵂ = invfourier(interp.position, kpoints)
    Āᴴ = map(zip(Aᵂ, U)) do (Aᵂₖ, Uₖ)
        Uₖ' * Aᵂₖ * Uₖ
    end

    # commutators [Dᴴα, Āᴴβ] - [Dᴴβ, Āᴴα]
    Dᴴ_Āᴴ = map(zip(Dᴴ, Āᴴ)) do (Dᴴₖ, Āᴴₖ)
        # each element: 3x1 vector * 1x3 vector -> 3x3 matrix,
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

    Ωᵂ = invfourier(interp.berry_curvature, kpoints)
    # gauge-covariant part of Berry curvature, Ω̄ᴴ
    Ω̄ᴴ = map(zip(Ωᵂ, U)) do (Ωᵂₖ, Uₖ)
        Uₖ' * Ωᵂₖ * Uₖ
    end

    Ωᴴ = Ω̄ᴴ .- Dᴴ_Āᴴ .- im .* Dᴴ_Dᴴ
    # just need band-diagonal elements
    return map(Ωᴴ) do Ωᴴₖ
        # force it to be real
        real(diag(Ωᴴₖ))
    end
end

"""Interpolate Berry curvature in Bloch gauge, sumed over bands with given
occupations, using WYSV06 Eq. 32."""
function (interp::BerryCurvatureInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}, ::WYSV06; kwargs...
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)

    eigenvalues, U, _, Dᴴ = compute_D_matrix(
        interp.hamiltonian, interp.hamiltonian_gradient, kpoints; kwargs...
    )
    occupations = [Int.(εₖ .<= interp.fermi_energy) for εₖ in eigenvalues]

    # WYSV Eq. 27
    # gauge-covariant part of k-space position operator
    Aᵂ = invfourier(interp.position, kpoints)
    Āᴴ = map(zip(Aᵂ, U)) do (Aᵂₖ, Uₖ)
        Uₖ' * Aᵂₖ * Uₖ
    end

    Ωᵂ = invfourier(interp.berry_curvature, kpoints)
    # gauge-covariant part of Berry curvature, Ω̄ᴴ
    Ω̄ᴴ = map(zip(Ωᵂ, U)) do (Ωᵂₖ, Uₖ)
        Uₖ' * Ωᵂₖ * Uₖ
    end

    return map(zip(Ω̄ᴴ, Āᴴ, Dᴴ, occupations)) do (Ω̄ᴴₖ, Āᴴₖ, Dᴴₖ, fₖ)
        # occupation -> matrix
        F = transpose(fₖ) .- fₖ
        DA = (F .* Dᴴₖ) * transpose.(Āᴴₖ)
        return real(
            sum(fₖ .* diag(Ω̄ᴴₖ)) +
            sum(diag(DA - transpose.(DA) + (im .* F .* Dᴴₖ) * transpose.(Dᴴₖ))),
        )
    end
end

"""Interpolate Berry curvature and sum over bands with given occupations,
using LVTS12 Eq. 51.

This is numerically more accurate than WYSV06, since it operates in Wannier
gauge and the spiky J matrix (= im * D matrix) is regulated by the occupation f:
f * J * g, where g = 1 - f.
"""
function (interp::BerryCurvatureInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}, ::LVTS12; kwargs...
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)

    # Wannier-gauge position operator, LVTS12 Eq. 28 and 78
    Aᵂ = invfourier(interp.position, kpoints)
    # Wannier-gauge Berry curvature, LVTS12 Eq. 31
    Ωᵂ = invfourier(interp.berry_curvature, kpoints)

    # Wannier-gauge, J⁻ = f*J*g, J⁺ = g*J*f, LVTS12 Eq. 76
    eigenvalues, U, _, Dᴴ = compute_D_matrix(
        interp.hamiltonian, interp.hamiltonian_gradient, kpoints; kwargs...
    )
    J, J⁻, J⁺ = compute_J_matrix(eigenvalues, U, Dᴴ, interp.fermi_energy)

    # Wannier-gauge occupation operator
    fᵂ = map(zip(eigenvalues, U)) do (εₖ, Uₖ)
        Uₖ * Diagonal(Int.(εₖ .<= interp.fermi_energy)) * Uₖ'
    end
    # LVTS12 Eq. 51
    return map(zip(fᵂ, Ωᵂ, Aᵂ, J, J⁻, J⁺)) do (fᵂₖ, Ωᵂₖ, Aᵂₖ, Jₖ, J⁻ₖ, J⁺ₖ)
        # note we use the cyclic property of trace for the 2nd term of the
        # RHS of LVTS12 Eq. 51, so that we can directly use J⁺ₖ
        real(tr(fᵂₖ * Ωᵂₖ)) -
        2 * imag(tr(Aᵂₖ * transpose.(J⁺ₖ) + J⁻ₖ * transpose.(Aᵂₖ + Jₖ)))
    end
end

function (interp::BerryCurvatureInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}; kwargs...
)
    return interp(kpoints, LVTS12(); kwargs...)
end

"""
    $(SIGNATURES)

Convert an axial vector to second-rank antisymmetric tensor.

`v_γ = ϵ_αβγ * v_αβ`

WYSV Eq. 5.

# Arguments
- `v`: axial vector, `v = (v_x, v_y, v_z)`
"""
function axialvector_to_antisymmetrictensor(v::AbstractArray)
    vx, vy, vz = v
    return MMat3([
        0 vz -vy
        -vz 0 vx
        vy -vx 0
    ])
end

"""
    $(SIGNATURES)

Convert a second-rank antisymmetric tensor to an axial vector.
"""
function antisymmetrictensor_to_axialvector(M::AbstractMatrix)
    return MVec3(M[2, 3], -M[1, 3], M[1, 2])
end
