export OrbitalMagnetizationInterpolator,
    TBHamiltonianPosition, TBPositionHamiltonianPosition

"""Construct a tight-binding <0| H (r - R) |R> operator (LVTS12 Eq. 83) operator."""
function TBHamiltonianPosition end

function TBHamiltonianPosition(Rspace::BareRspace, operator::AbstractVector)
    @assert !isempty(operator) "empty operator"
    @assert !isempty(operator[1]) "empty operator"
    @assert operator[1] isa AbstractMatrix "operator must be a matrix"
    v = operator[1][1, 1]
    @assert v isa AbstractVector && length(v) == 3 "each element must be 3-vector"
    T = real(eltype(v))
    M = Matrix{MVec3{Complex{T}}}
    return TBOperator{M}("HamiltonianPosition", Rspace, operator)
end

"""
    $(SIGNATURES)

# Arguments
- `Rspace`: should be a [`WignerSeitzRspace`](@ref) or [`MDRSRspace`](@ref), and
    it will be simplified to a [`BareRspace`](@ref) by [`simplify`](@ref).
- `overlaps`, `eigenvalues`: should be Hamiltonian-gauge matrices, `n_bands`-sized
- `gauges`: from Hamiltonian gauge to Wannier gauge, `n_bands * n_wannier`
"""
function TBHamiltonianPosition(
    Rspace::Union{WignerSeitzRspace,MDRSRspace},
    kstencil::KspaceStencil,
    overlaps::AbstractVector,
    eigenvalues::AbstractVector,
    gauges::AbstractVector,
)
    # LVTS12 Eq. 91
    # H_kb = U'_k H_k M_kb U_kb
    kpb_k = kstencil.kpb_k
    kpb_G = kstencil.kpb_G
    Hkb = compute_hamiltonian_times_position_kspace(overlaps, kpb_k, kpb_G, eigenvalues)
    # Hkb transforms like overlap matrices
    Hkbᵂ = transform_gauge(Hkb, kpb_k, gauges)

    # sum over bvectors, RHS of LVTS12 Eq. 83:
    # <0| H (r - R) |R> = 1/Nₖ * i * ∑_kb exp(-ikR) w_b b H_kb
    nkpts = n_kpoints(kstencil)
    kpoints = kstencil.kpoints
    recip_lattice = reciprocal_lattice(kstencil)
    wb = kstencil.bweights
    nbvecs = n_bvectors(kstencil)
    T = eltype(Hkbᵂ[1][1])
    nwann = size(Hkbᵂ[1][1], 1)

    Hr_k = map(1:nkpts) do ik
        Hrₖ = zeros(Vec3{T}, nwann, nwann)
        for ib in 1:nbvecs
            ik2 = kpb_k[ik][ib]
            G = kpb_G[ik][ib]
            # b isa Vec3, along x, y, z directions
            b = recip_lattice * (kpoints[ik2] + G - kpoints[ik])
            Hrₖ .+= im * wb[ib] .* Ref(b) .* Hkbᵂ[ik][ib]
        end
        return Hrₖ
    end
    # Hr_k indexed as Hr_k[ik][m, n][α]

    Hr_R = fourier(kpoints, Hr_k, Rspace)
    bare_Rspace, bare_Hr_R = simplify(Rspace, Hr_R)
    return TBHamiltonianPosition(bare_Rspace, bare_Hr_R)
end

function TBHamiltonianPosition(
    Rspace::Union{WignerSeitzRspace,MDRSRspace},
    model::Model,
    gauges::AbstractVector=model.gauges,
)
    return TBHamiltonianPosition(
        Rspace, model.kstencil, model.overlaps, model.eigenvalues, gauges
    )
end

function TBHamiltonianPosition(model::Model, gauges::AbstractVector=model.gauges; kwargs...)
    Rspace = generate_Rspace(model; kwargs...)
    return TBHamiltonianPosition(Rspace, model, gauges)
end

"""Construct a tight-binding <0| H (r - R) |R> operator (LVTS12 Eq. 83) operator."""
function TBPositionHamiltonianPosition end

function TBPositionHamiltonianPosition(Rspace::BareRspace, operator::AbstractVector)
    @assert !isempty(operator) "empty operator"
    @assert !isempty(operator[1]) "empty operator"
    @assert operator[1] isa AbstractMatrix "operator must be a matrix"
    v = operator[1][1, 1]
    @assert v isa AbstractMatrix && size(v) == (3, 3) "each element must be 3*3 matrix"
    T = real(eltype(v))
    M = Matrix{MMat3{Complex{T}}}
    return TBOperator{M}("PositionHamiltonianPosition", Rspace, operator)
end

"""
    $(SIGNATURES)

# Arguments
- `Rspace`: should be a [`WignerSeitzRspace`](@ref) or [`MDRSRspace`](@ref), and
    it will be simplified to a [`BareRspace`](@ref) by [`simplify`](@ref).
- `overlaps`, `eigenvalues`: should be Hamiltonian-gauge matrices, `n_bands`-sized
- `gauges`: from Hamiltonian gauge to Wannier gauge, `n_bands * n_wannier`
"""
function TBPositionHamiltonianPosition(
    Rspace::Union{WignerSeitzRspace,MDRSRspace},
    kstencil::KspaceStencil,
    uHu::AbstractVector,
    gauges::AbstractVector;
    force_hermiticity=default_w90_berry_duHdu_force_hermiticity(),
)
    nkpts = n_kpoints(kstencil)
    nbvecs = n_bvectors(kstencil)
    # LVTS12 Eq. 93
    # H_k,b1,b2ᵂ = U_k,b1' H_k,b1,b2 U_k,b2
    uHuᵂ = map(1:nkpts) do ik
        map(CartesianIndices((nbvecs, nbvecs))) do idx
            ib1, ib2 = idx.I
            ikb1 = kstencil.kpb_k[ik][ib1]
            ikb2 = kstencil.kpb_k[ik][ib2]
            return gauges[ikb1]' * uHu[ik][ib1, ib2] * gauges[ikb2]
        end
    end

    # sum over bvectors, RHS of LVTS12 Eq. 84:
    # <0| r H (r - R) |R> = 1/Nₖ * ∑_k,b1,b2 exp(-ikR) w_b1 b1 w_b2 b2 H_k,b1,b2
    kpoints = kstencil.kpoints
    recip_lattice = reciprocal_lattice(kstencil)
    wb = kstencil.bweights
    T = eltype(uHuᵂ[1][1])
    nwann = size(uHuᵂ[1][1], 1)

    rHr_k = map(1:nkpts) do ik
        rHrₖ = zeros(MMat3{T}, nwann, nwann)
        for ib2 in 1:nbvecs
            for ib1 in 1:nbvecs
                ikb1 = kstencil.kpb_k[ik][ib1]
                ikb2 = kstencil.kpb_k[ik][ib2]
                G1 = kstencil.kpb_G[ik][ib1]
                G2 = kstencil.kpb_G[ik][ib2]
                # b isa Vec3, along x, y, z directions
                b1 = recip_lattice * (kpoints[ikb1] + G1 - kpoints[ik])
                b2 = recip_lattice * (kpoints[ikb2] + G2 - kpoints[ik])
                rHrₖ .+= (wb[ib1] * wb[ib2]) .* Ref(b1 * b2') .* uHuᵂ[ik][ib1, ib2]
            end
        end
        return rHrₖ
    end
    # rHr_k indexed as rHr_k[ik][m, n][α, β]
    force_hermiticity && hermitize_duHdu!(rHr_k)

    rHr_R = fourier(kpoints, rHr_k, Rspace)
    bare_Rspace, bare_Hr_R = simplify(Rspace, rHr_R)
    return TBPositionHamiltonianPosition(bare_Rspace, bare_Hr_R)
end

"""
Make `< ∂α uₘₖ | Hₖ | ∂β uₙₖ >` Hermitian.

In theory, `rHrₖ[m, n][α, β] = < ∂α uₘₖ | Hₖ | ∂β uₙₖ >` should satisfy
Hermiticity of `< ∂α uₘₖ | Hₖ | ∂β uₙₖ >† = < ∂β uₙₖ | Hₖ | ∂α uₘₖ >`;
Here we enforce this to be safe.
"""
function hermitize_duHdu!(duHdu::AbstractVector)
    nkpts = length(duHdu)
    nwann = size(duHdu[1], 1)
    for ik in 1:nkpts
        for n in 1:nwann
            for m in 1:nwann
                for α in 1:3
                    for β in 1:α
                        duHdu[ik][m, n][α, β] = conj(duHdu[ik][n, m][β, α])
                    end
                end
            end
        end
    end
    return nothing
end

function TBPositionHamiltonianPosition(
    Rspace::Union{WignerSeitzRspace,MDRSRspace},
    model::Model,
    uHu::AbstractVector,
    gauges::AbstractVector=model.gauges,
)
    return TBPositionHamiltonianPosition(Rspace, model.kstencil, uHu, gauges)
end

function TBPositionHamiltonianPosition(
    model::Model, uHu::AbstractVector, gauges::AbstractVector=model.gauges
)
    Rspace = generate_Rspace(model)
    return TBPositionHamiltonianPosition(Rspace, model, uHu, gauges)
end

"""
    $(TYPEDEF)

A struct for interpolating orbital magnetization on given kpoints.

# Fields
$(FIELDS)
"""
struct OrbitalMagnetizationInterpolator <: AbstractTBInterpolator
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

    """R-space <0| H (r - R) |R> operator, LVTS12 Eq. 83."""
    hamiltonian_position::TBOperator

    """R-space <0| r H (r - R) |R> operator, LVTS12 Eq. 83."""
    position_hamiltonian_position::TBOperator

    """Fermi energy in eV"""
    fermi_energy::Float64
end

function OrbitalMagnetizationInterpolator(
    hamiltonian::TBOperator,
    position::TBOperator,
    hamiltonian_position::TBOperator,
    position_hamiltonian_position::TBOperator,
    fermi_energy::Real,
)
    return OrbitalMagnetizationInterpolator(
        hamiltonian,
        TBHamiltonianGradient(hamiltonian),
        position,
        TBBerryCurvature(position),
        hamiltonian_position,
        position_hamiltonian_position,
        fermi_energy,
    )
end

"""Interpolate orbital magnetization, LVTS Eq. 71 and 72."""
function (interp::OrbitalMagnetizationInterpolator)(
    kpoints::AbstractVector{<:AbstractVector};
    force_hermiticity=default_w90_berry_duHdu_force_hermiticity(),
    kwargs...,
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)

    # Wannier-gauge position operator, LVTS12 Eq. 28 and 78
    Aᵂ = invfourier(interp.position, kpoints)
    # Wannier-gauge Berry curvature, LVTS12 Eq. 31
    Ωᵂ = invfourier(interp.berry_curvature, kpoints)
    # Wannier-gauge Hamiltonian * position, LVTS12 Eq. 33 and 79
    Bᵂ = invfourier(interp.hamiltonian_position, kpoints)
    # Wannier-gauge position * Hamiltonian * position, LVTS12 Eq. 35 and 79
    Cᵂ = invfourier(interp.position_hamiltonian_position, kpoints)
    force_hermiticity && hermitize_duHdu!(Cᵂ)
    # Wannier-gauge Λ = i C - i C†, LVTS12 Eq. 37
    # Hermitian conjugate on the band index, the inner α, β indices are not transposed
    Λᵂ = im * (Cᵂ - permutedims.(conj(Cᵂ)))
    # Wannier-gauge Hamiltonian
    Hᵂ = invfourier(interp.hamiltonian, kpoints)
    # Rα * < m0 | H | nR >
    RHᵂ = invfourier(interp.hamiltonian_gradient, kpoints)
    eigenvalues, U, _, Dᴴ = compute_D_matrix(Hᵂ, RHᵂ, kpoints; kwargs...)
    # Wannier-gauge, J⁻ = f*J*g, J⁺ = g*J*f, LVTS12 Eq. 76
    εF = interp.fermi_energy
    _, J⁻, J⁺ = compute_J_matrix(eigenvalues, U, Dᴴ, εF)
    # Wannier-gauge occupation operator
    fᵂ = map(zip(eigenvalues, U)) do (εₖ, Uₖ)
        Uₖ * Diagonal(Int.(εₖ .<= εF)) * Uₖ'
    end

    # The total orbital magnetization can be decomposed into two parts:
    # - Local circulation (LC), LVTS Eq. 4
    #   M̃ᴸ = - e/(2ħ) ∫dk -2 Im(G - εF * F)
    # - Itinerant circulation (IC), LVTS Eq. 5
    #   M̃ᴵ = - e/(2ħ) ∫dk -2 Im(H - εF * F)
    return map(zip(fᵂ, Hᵂ, Ωᵂ, Aᵂ, Bᵂ, Λᵂ, J⁻, J⁺)) do (f, Hₖ, Ωₖ, Aₖ, Bₖ, Λₖ, J⁻ₖ, J⁺ₖ)
        H⁰ = f * Hₖ * f
        Ω⁰ = f * Ωₖ * f
        A⁰ = f * Aₖ * f
        Λ⁰ = f * Λₖ * f
        g = I - f
        H¹ = g * Hₖ * g
        A⁻ = f * Aₖ * g
        A⁺ = g * Aₖ * f
        B⁺ = g * Bₖ * f
        # Berry curvature, F := -2 Im(Fαβ), LVTS12 Eq. 51
        AJJAJ = A⁻ * transpose.(J⁺ₖ) + J⁻ₖ * transpose.(A⁺ + J⁺ₖ)
        F = real(tr(Ω⁰)) - 2 * imag(tr(AJJAJ))
        # IC = H - εF * F, where H := -2 Im(Hαβ), LVTS12 Eq. 56, also Eq. 71
        HAA = 2 * imag(tr(H⁰ * A⁰ * transpose.(A⁰)))
        H = real(tr(H⁰ * Ω⁰)) + HAA - 2 * imag(tr(H⁰ * AJJAJ))
        # LC = G - εF * F, where G := -2 Im(Gαβ), LVTS12 Eq. 66, also Eq. 72
        # note that the α ⟷ β exchange in LVTS12 Eq. 66 (and 72) only transpose
        # the α, β index, while the band indices m, n are unchanged
        JB = J⁻ₖ * transpose.(B⁺)
        JBJB = JB - transpose.(JB)
        G = real(tr(Λ⁰)) - HAA - 2 * imag(tr(JBJB + J⁻ₖ * H¹ * transpose.(J⁺ₖ)))
        return H + G - 2 * εF * F
    end
end
