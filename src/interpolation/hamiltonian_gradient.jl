"""
    $(TYPEDEF)

A struct for interpolating 1st-order derivative of Hamiltonian on given kpoints.

YWVS Eq. 26.

# Fields
$(FIELDS)
"""
struct HamiltonianGradientInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian"""
    hamiltonian::TBOperator
end

"""
Compute the derivative of the Hamiltonian with respect to three Cartesian
directions.

YWVS Eq. 26.
"""
function (interp::HamiltonianGradientInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}; kwargs...
)
    # since `KPathInterpolant` is also a `AbstractVector{<:AbstractVector}`,
    # we call the `get_kpoints` so that this function also works for `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    eigvals, _, dH, D_matrices = compute_D_matrix(interp.hamiltonian, kpoints; kwargs...)

    # dH is the gauge-covariant part, now compute non-diagonal part, to build the full dHᴴ
    # Here we allocate a buffer for permutedims!, which is faster than permutedims
    # since the latter allocates memory
    HDperm = similar(D_matrices[1])
    dHᴴ = map(zip(eigvals, dH, D_matrices)) do (εₖ, dHₖ, Dₖ)
        HD = Diagonal(εₖ) * Dₖ
        # Note we need permutedims because we only want to swap the WF index,
        # while adjoint() will recursively tranpose the inner MVec3.
        # The conj() instead should apply recursively to the inner MVec3.
        permutedims!(HDperm, HD, [2, 1])
        dHₖ + HD + conj(HDperm)
    end
    return dHᴴ
end

"""
    $(TYPEDEF)

A struct for interpolating velocity on given kpoints.

Velocity is the diagonal part of the derivative of Hamiltonian in Bloch gauge,
a simpler (and a little bit faster) version of
[`HamiltonianGradientInterpolator`](@ref).

YWVS Eq. 27.

# Fields
$(FIELDS)
"""
struct VelocityInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian"""
    hamiltonian::TBOperator
end

abstract type AbstractVelocityAlgorithm end

"""Compute velocity using inverse Fourier transform of ``\\mathbf{R} H`` operator."""
struct FourierSpaceVelocity <: AbstractVelocityAlgorithm end

"""Compute velocity using finite difference of Wannier-interpolated eigenvalues."""
struct FiniteDifferenceVelocity <: AbstractVelocityAlgorithm end

"""
Compute velocity (in Bloch gauge) along three Cartesian directions.

# Return
- `eigenvalues`: energy eigenvalues
- `velocities`: velocity along three Cartesian directions, in unit `hbar * m / s`

!!! warning

    `Wannier90` by default set `use_degen_pert = false`.
    In 3D, and for N degenerate states, the velocity is a tensor
    of size N * N * 3, where 3 is for three Cartesian directions.
    Thus I cannot simultaneously diagonalize the tensor for all 3 directions.
    This means I can only use perturbation treatment for one of the directions,
    and only in that direction the velocity matrix is diagonal.
    So for degenerate states, the velocity is not well defined, and the results
    are meaningless, instead one should use the full velocity matrix which
    also include non-diagonal part, see [`HamiltonianGradientInterpolator`](@ref).
"""
function (interp::VelocityInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}, ::FourierSpaceVelocity; kwargs...
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    _, _, dH, _ = compute_D_matrix(interp.hamiltonian, kpoints; kwargs...)
    # velocity is the diagonal part, not affected by the D matrices
    return map(dH) do dHₖ
        # YWVS Eq.27
        real(diag(dHₖ))
    end
end

"""
Compute the velocity using finite differences of 2nd order.

PRB 93, 205147 (2016)  Eq. 80.

# Keyword arguments
- `dk`: the finite difference spacing, Å⁻¹ unit, default to `1e-3`
"""
function (interp::VelocityInterpolator)(
    kpoints::AbstractVector{<:AbstractVector},
    ::FiniteDifferenceVelocity;
    dk::Real=1e-3,
    kwargs...,
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    nwann = n_wannier(interp.hamiltonian)
    nkpts = length(kpoints)

    # the computed velocity along 3 Cartesian directions
    @assert n_Rvectors(interp.hamiltonian) > 0 "empty Hamiltonian"
    T = eltype(interp.hamiltonian[1])
    RT = real(T)
    velocity = [zeros_eigenvalues(RT, nkpts, nwann) for _ in 1:3]
    # the 6 interpolated kpoints, each element is a Cartesian coordinates
    Δk⁺ = Vec3{T}[[dk, 0, 0], [0, dk, 0], [0, 0, dk]]
    Δk⁻ = -Δk⁺

    # the interpolated Hamiltonain
    H_k⁺ = zeros_gauge(T, nkpts, nwann)
    H_k⁻ = zeros_gauge(T, nkpts, nwann)
    # the interpolated eigenvalues
    eigvals⁺ = zeros_eigenvalues(RT, nkpts, nwann)
    eigvals⁻ = zeros_eigenvalues(RT, nkpts, nwann)
    # buffer for eigenvectors
    eigvecs⁺ = zeros_gauge(T, nkpts, nwann)
    eigvecs⁻ = zeros_gauge(T, nkpts, nwann)
    # to fractional
    inv_recip_latt = inv(reciprocal_lattice(real_lattice(interp.hamiltonian)))
    Δk⁺_frac = map(v -> inv_recip_latt * v, Δk⁺)
    Δk⁻_frac = map(v -> inv_recip_latt * v, Δk⁻)

    # three Cartesian directions
    for idir in 1:3
        invfourier!(H_k⁺, interp.hamiltonian, kpoints .+ Ref(Δk⁺_frac[idir]))
        invfourier!(H_k⁻, interp.hamiltonian, kpoints .+ Ref(Δk⁻_frac[idir]))

        eigen!(eigvals⁺, eigvecs⁺, H_k⁺)
        eigen!(eigvals⁻, eigvecs⁻, H_k⁻)

        velocity[idir] .= (eigvals⁺ .- eigvals⁻) ./ (2dk)
    end

    # convert to nested MVec3, result can be accessed as `velocity[ik][iw]`
    return map(1:nkpts) do ik
        map(1:nwann) do iw
            MVec3{RT}(velocity[1][ik][iw], velocity[2][ik][iw], velocity[3][ik][iw])
        end
    end
end

function (interp::VelocityInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}; kwargs...
)
    return interp(kpoints, FourierSpaceVelocity(); kwargs...)
end
