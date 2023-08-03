"""
    $(TYPEDEF)

A struct for interpolating 2nd-order derivative of Hamiltonian on given kpoints.

YWVS Eq. 28.

# Fields
$(FIELDS)
"""
struct HamiltonianHessianInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian"""
    hamiltonian::TBOperator
end

"""
    $(SIGNATURES)

Compute the second derivative of the Hamiltonian with respect to three
Cartesian directions.

YWVS Eq. 28.
"""
function (interp::HamiltonianHessianInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}; kwargs...
)
    @assert n_Rvectors(interp.hamiltonian) > 0 "empty Hamiltonian"
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)

    _, U, dH, D_matrices = compute_D_matrix(interp.hamiltonian, kpoints; kwargs...)

    # Objective: 2nd order derivative d²H := ∂ᵦ Vα = ∂²H / ∂kα ∂kβ
    # The following comments use notations:
    # - α, β ∈ {x, y, z} for three Cartesian directions
    # - m, n ∈ 1:nwann for Wannier functions

    # to Cartesian, in Å
    R_cart = map(v -> interp.hamiltonian.lattice * v, interp.hamiltonian.Rvectors)
    # alias for convenience
    H_R = interp.hamiltonian.operator

    # 1. Construct Rspace operator: -Rα * Rβ * H(R), YWVS Eq. 30
    # 1.1. -Rα * Rβ, index by RR[iR][α, β]
    RR = -R_cart .* transpose.(R_cart)
    # 1.2. -Rα * Rβ * H(R), index by RRH_R[iR][m, n][α, β]
    RRH_R = map(zip(RR, H_R)) do (rr, h)
        # make rr broadcastable to h, the result can be indexed by [m, n][α, β]
        Ref(rr) .* h
    end
    # 2. inv Fourier to k-space, indexed by d²Hᵂ[ik][m, n][α, β]
    d²Hᵂ = invfourier(interp.hamiltonian.Rspace, RRH_R, kpoints)
    # 3. to Bloch gauge, indexed by d²Hᴴ[ik][m, n][α, β]
    # 3.1. the covariant part of d²H, i.e., the
    #      ``\bar{H}_{\alpha\beta}^{(H)}`` in YWVS Eq. 28
    d²Hᴴ = adjoint.(U) .* d²Hᵂ .* U
    # 3.2. add remaing part for gauge transformation, YWVS Eq. 28
    HD = map(zip(dH, D_matrices)) do (dHₖ, Dₖ)
        # dHₖ and Dₖ are both Matrix{MVec3}
        # transpose the inner MVec3 of Dₖ so that the returned result is
        # indexed by [m, n][α, β]
        dHₖ * transpose.(Dₖ)
    end
    # note the conjugate transpose is on the [m, n] index, not the α,β index
    # need to use permutedims to swap the WF indinces [m, n], because transpose()
    # is recursive
    d²Hᴴ .+= HD .+ conj(permutedims.(HD))

    # the final result is indexed by d²Hᴴ[ik][m, n][α, β]
    return d²Hᴴ
end

"""
    $(TYPEDEF)

A struct for interpolating effective mass on given kpoints.

i.e., the diagonal part of 2nd-order derivative of Hamiltonian.

YWVS Eq. 28.

# Fields
$(FIELDS)
"""
struct EffectiveMassInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian"""
    hamiltonian::TBOperator
end

abstract type AbstractEffectiveMassAlgorithm end

"""Compute effective mass using inverse Fourier transform of ``\\mathbf{R} H`` operator."""
struct AnalyticEffectiveMass <: AbstractEffectiveMassAlgorithm end

"""Compute effective mass using finite difference of Wannier-interpolated eigenvalues."""
struct FiniteDifferenceEffectiveMass <: AbstractEffectiveMassAlgorithm end

@inline function (interp::EffectiveMassInterpolator)(
    kpi::KPathInterpolant, args...; kwargs...
)
    kpoints = get_kpoints(kpi)
    return interp(kpoints, args...; kwargs...)
end

"""
    $(SIGNATURES)

Compute the inverse of effective mass using Wannier interpolation.

YWVS Eq.28
"""
function (interp::EffectiveMassInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}, ::AnalyticEffectiveMass; kwargs...
)
    hessian_interp = HamiltonianHessianInterpolator(interp.hamiltonian)
    d²Hᴴ = hessian_interp(kpoints; kwargs...)

    return map(d²Hᴴ) do d²Hₖ
        # YWVS Eq.28
        real(diag(d²Hₖ))
    end
end

"""
    $(SIGNATURES)

Compute the inverse of effective mass using finite differences of 2nd order.

Apply twice PRB 93, 205147 (2016)  Eq. 80.

# Keyword arguments
- `dk`: the finite difference spacing, Å⁻¹ unit, default to `1e-3`
"""
function (interp::EffectiveMassInterpolator)(
    kpoints::AbstractVector{<:AbstractVector},
    ::FiniteDifferenceEffectiveMass;
    dk=1e-3,
    kwargs...,
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    nwann = n_wannier(interp.hamiltonian)
    nkpts = length(kpoints)
    @assert n_Rvectors(interp.hamiltonian) > 0 "empty Hamiltonian"
    T = real(eltype(interp.hamiltonian[1]))

    # the 1 + 6 + 12 interpolated kpoints, each column is a fractional coordinates
    # note the original kpoint is included
    # 6 points at k ± dk, for computing dε/dk at k ± dk/2, along the same direction as dk
    # 12 points at k ± dk/2, for computing dε/dk at k ± dk/2, but along the other 2 directions
    Δk = Vec3[
        [0, 0, 0],
        [-dk, 0, 0],
        [dk, 0, 0],
        [0, -dk, 0],
        [0, dk, 0],
        [0, 0, -dk],
        [0, 0, dk],
        [-dk / 2, -dk / 2, 0],
        [-dk / 2, dk / 2, 0],
        [-dk / 2, 0, -dk / 2],
        [-dk / 2, 0, dk / 2],
        [dk / 2, -dk / 2, 0],
        [dk / 2, dk / 2, 0],
        [dk / 2, 0, -dk / 2],
        [dk / 2, 0, dk / 2],
        [0, -dk / 2, -dk / 2],
        [0, -dk / 2, dk / 2],
        [0, dk / 2, -dk / 2],
        [0, dk / 2, dk / 2],
    ]
    # to fractional
    inv_recip_latt = inv(reciprocal_lattice(real_lattice(interp.hamiltonian)))
    Δk_frac = map(k -> inv_recip_latt * k, Δk)

    # alias
    H_R = interp.hamiltonian.operator
    # the interpolated Hamiltonain, indexed by H_k[iδk][ik][m, n], where
    # - iδk ∈ 1:length(Δk_frac) is the index of Δk_frac
    # - ik for kpoints
    # - m, n for Wannier functions
    H_k = map(Δk_frac) do δk
        invfourier(interp.hamiltonian.Rspace, H_R, kpoints .+ Ref(δk))
    end
    # the interpolated eigenvalues, indexed by ε[iδk][ik][n]
    ε = map(H_k) do H
        eigen(H)[1]
    end
    # dε/dk at 6 points: k ± dk/2, along 3 Cartesian directions
    # indexed by dε[iδk][α][ik][n], where
    # - iδk ∈ 1:6 for the first 6 kpoints of Δk_frac
    # - α ∈ 1:3 for Cartesian directions
    dε = [[zeros_eigenvalues(T, nkpts, nwann) for _ in 1:3] for _ in 1:6]
    # dε/dk at k ± dk/2, along the same direction as dk
    dε[1][1] .= (ε[2] .- ε[1]) ./ -dk
    dε[2][1] .= (ε[3] .- ε[1]) ./ dk
    dε[3][2] .= (ε[4] .- ε[1]) ./ -dk
    dε[4][2] .= (ε[5] .- ε[1]) ./ dk
    dε[5][3] .= (ε[6] .- ε[1]) ./ -dk
    dε[6][3] .= (ε[7] .- ε[1]) ./ dk

    # dε/dk at k ± dk/2, but along the other 2 directions
    dε[1][2] .= (ε[9] .- ε[8]) ./ dk
    dε[1][3] .= (ε[11] .- ε[10]) ./ dk
    dε[2][2] .= (ε[13] .- ε[12]) ./ dk
    dε[2][3] .= (ε[15] .- ε[14]) ./ dk
    dε[3][1] .= (ε[12] .- ε[8]) ./ dk
    dε[3][3] .= (ε[17] .- ε[16]) ./ dk
    dε[4][1] .= (ε[13] .- ε[9]) ./ dk
    dε[4][3] .= (ε[19] .- ε[18]) ./ dk
    dε[5][1] .= (ε[14] .- ε[10]) ./ dk
    dε[5][2] .= (ε[18] .- ε[16]) ./ dk
    dε[6][1] .= (ε[15] .- ε[11]) ./ dk
    dε[6][2] .= (ε[19] .- ε[17]) ./ dk

    # d²ε/dk² at each kpoints[ik], along 3 Cartesian directions
    # indexed by d²ε[α][β][ik][n], where α, β ∈ 1:3 for Cartesian directions
    d²ε = [[zeros_eigenvalues(T, nkpts, nwann) for _ in 1:3] for _ in 1:3]
    # d²ε/dk² at k
    d²ε[1][1] .= (dε[2][1] .- dε[1][1]) ./ dk
    d²ε[1][2] .= (dε[4][1] .- dε[3][1]) ./ dk
    d²ε[1][3] .= (dε[6][1] .- dε[5][1]) ./ dk
    d²ε[2][1] .= (dε[2][2] .- dε[1][2]) ./ dk
    d²ε[2][2] .= (dε[4][2] .- dε[3][2]) ./ dk
    d²ε[2][3] .= (dε[6][2] .- dε[5][2]) ./ dk
    d²ε[3][1] .= (dε[2][3] .- dε[1][3]) ./ dk
    d²ε[3][2] .= (dε[4][3] .- dε[3][3]) ./ dk
    d²ε[3][3] .= (dε[6][3] .- dε[5][3]) ./ dk

    # nest α, β
    return map(1:nkpts) do ik
        map(1:nwann) do n
            MMat3{T}(
                [
                    d²ε[1][1][ik][n] d²ε[1][2][ik][n] d²ε[1][3][ik][n]
                    d²ε[2][1][ik][n] d²ε[2][2][ik][n] d²ε[2][3][ik][n]
                    d²ε[3][1][ik][n] d²ε[3][2][ik][n] d²ε[3][3][ik][n]
                ],
            )
        end
    end
end
