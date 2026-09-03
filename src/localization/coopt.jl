export SpinModel

"""
Model for spin polarized system with constraint.

Traditionally, we run two independent Wannierizations for spin up and spin down.
Here we add a constraint to maximally overlap the spin-up and spin-down WFs,
so that they map one-by-one to each other.

The constructor enforces that `up` and `dn` describe the same real-space /
reciprocal-space system (lattice, atoms, kstencil); the `overlaps_updn` field stores the
Bloch-basis ↑↓ overlap `< u_nk^↑ | u_mk^↓ >` read from the coupling amn file.
"""
struct SpinModel{T <: Real}
    up::Model{T}
    dn::Model{T}

    "< u_nk^↑ | u_mk^↓ >, size: (n_bands, n_bands, n_kpts)"
    overlaps_updn::Array{Complex{T}, 3}

    function SpinModel{T}(up::Model{T}, dn::Model{T}, M::AbstractArray{<:Complex, 3}) where {T <: Real}
        isapprox(up.lattice, dn.lattice) || error("SpinModel: up/dn lattice mismatch")
        (
            length(up.atom_positions) == length(dn.atom_positions) &&
                all(isapprox(a, b) for (a, b) in zip(up.atom_positions, dn.atom_positions))
        ) ||
            error("SpinModel: up/dn atom_positions mismatch")
        up.atom_labels == dn.atom_labels || error("SpinModel: up/dn atom_labels mismatch")
        isapprox(up.kstencil, dn.kstencil) || error("SpinModel: up/dn kstencil mismatch")
        size(M, 3) == n_kpoints(up) || error("SpinModel: M has wrong number of kpoints")
        size(M, 1) == size(M, 2) == n_bands(up) || error("SpinModel: M has wrong n_bands")
        return new{T}(up, dn, Array{Complex{T}, 3}(M))
    end
end

function SpinModel(up::Model{T}, dn::Model{T}, M::AbstractArray{<:Complex, 3}) where {T <: Real}
    return SpinModel{T}(up, dn, M)
end

function Base.show(io::IO, ::MIME"text/plain", model::SpinModel)
    println(io, "spin up:")
    show(io, model.up)
    println(io, "\n")

    println(io, "spin down:")
    return show(io, model.dn)
end

# Accessors delegate to `up` — constructor invariants guarantee they match `dn`.
n_atoms(m::SpinModel) = n_atoms(m.up)
n_kpoints(m::SpinModel) = n_kpoints(m.up)
n_bvectors(m::SpinModel) = n_bvectors(m.up)
n_bands(m::SpinModel) = n_bands(m.up)
n_wannier(m::SpinModel) = n_wannier(m.up)
CrystalBase.real_lattice(m::SpinModel) = real_lattice(m.up)
CrystalBase.reciprocal_lattice(m::SpinModel) = reciprocal_lattice(m.up)
kpoints(m::SpinModel) = kpoints(m.up)
kgrid_size(m::SpinModel) = kgrid_size(m.up)
kpb_k(m::SpinModel) = kpb_k(m.up)
kpb_G(m::SpinModel) = kpb_G(m.up)
bweights(m::SpinModel) = bweights(m.up)

struct SpreadMag{T <: Real, S <: AbstractSpread} <: AbstractSpread
    # up spread
    up::S

    # dn spread
    dn::S

    # unit Å²
    Ωupdn::T

    # Total spread Ωt = up.Ω + dn.Ω + λ * Ω↑↓
    Ωt::T

    # overlap matrix between up and down WFs, unit Å², size = (n_wann, n_wann)
    M::Matrix{T}

    # λ
    λ::T
end

function spread(model::SpinModel, Uup, Udn, λ::Real)
    up = spread(model.up, Uup)
    dn = spread(model.dn, Udn)
    M = overlap_updn(model, Uup, Udn)
    Ωupdn = omega_updn(M)
    Ωt = up.Ω + dn.Ω + λ * Ωupdn
    return SpreadMag(up, dn, Ωupdn, Ωt, M, λ)
end

function spread(model::SpinModel, λ::Real)
    return spread(model, model.up.gauges, model.dn.gauges, λ)
end

function Base.show(io::IO, M::MIME"text/plain", Ω::SpreadMag)
    @info "spin up:"
    show(io, M, Ω.up)
    println(io, "\n")

    @info "spin down:"
    show(io, M, Ω.dn)
    println(io, "\n")

    n_wann = size(Ω.M, 1)
    @info "overlap between up and down WFs:"
    @printf(io, "  WF     <↑|↓>/Å²\n")
    for i in 1:n_wann
        @printf(io, "%4d %11.5f\n", i, Ω.M[i, i])
    end
    println(io, "")

    @info "Sum spread: Ωt = Ω↑ + Ω↓ + λ * Ω↑↓"
    @printf(io, "   Ω↑  = %11.5f\n", Ω.up.Ω)
    @printf(io, "   Ω↓  = %11.5f\n", Ω.dn.Ω)
    @printf(io, "   Ω↑↓ = %11.5f\n", Ω.Ωupdn)
    return @printf(io, "   Ωt  = %11.5f\n", Ω.Ωt)
end

"""
    overlap_updn(M, Uup, Udn)

Compute the overlap between up and down WFs.

Actually N - Ω↑↓, according to QPPM Eq. 8, where N = n_wann.

# Arguments
- `M`: the `SpinModel.overlaps_updn` matrices, `(n_bands, n_bands, n_kpts)`
- `Uup`: the up gauge matrices, `(n_bands, n_wann, n_kpts)`
- `Udn`: the down gauge matrices, `(n_bands, n_wann, n_kpts)`
"""
function overlap_updn(
        M::AbstractArray{<:Complex, 3},
        Uup::AbstractArray{<:Complex, 3},
        Udn::AbstractArray{<:Complex, 3},
    )
    n_bands, n_wann, n_kpts = size(Uup)
    Mᵂ = zeros(eltype(M), n_wann, n_wann)

    for ik in 1:n_kpts
        Mᵂ .+= view(Uup, :, :, ik)' * view(M, :, :, ik) * view(Udn, :, :, ik)
    end

    return abs2.(Mᵂ) ./ n_kpts^2
end

function overlap_updn(model::SpinModel, Uup, Udn)
    return overlap_updn(model.overlaps_updn, Uup, Udn)
end

function overlap_updn(model::SpinModel)
    return overlap_updn(model, model.up.gauges, model.dn.gauges)
end

"""
    omega_updn(M)

Compute QPPM Eq. 8.

# Arguments
- `M`: the overlap matrix between up and down WFs, size: (n_wann, n_wann),
    should be the matrix returned from [`overlap_updn`](@ref).
"""
function omega_updn(M::AbstractMatrix)
    n_wann = size(M, 1)
    # I am using minus sign here because the optimizer is minimizing total spread,
    # thus maximizing the ↑↓ overlap.
    return n_wann - sum(diag(M))
end

function omega_updn(model::SpinModel, Uup, Udn)
    return omega_updn(overlap_updn(model, Uup, Udn))
end

function omega_updn(model::SpinModel)
    return omega_updn(overlap_updn(model))
end

# Overload accepting M directly (for solver workspace use)
function omega_updn(M::AbstractArray{<:Complex, 3}, Uup::AbstractArray{<:Complex, 3}, Udn::AbstractArray{<:Complex, 3})
    return omega_updn(overlap_updn(M, Uup, Udn))
end

@doc raw"""
    overlap_updn_grad(model::SpinModel, Uup, Udn)

Compute gradients of [`overlap_updn`](@ref overlap_updn).

``\frac{d \Omega}{d U^{\uparrow}}`` and ``\frac{d \Omega}{d U^{\downarrow}}``.

TODO: this is actually the gradient of Tr[overlap_updn]

# Arguments
- `M`: the `SpinModel.overlaps_updn` matrices, `(n_bands, n_bands, n_kpts)`
- `Uup`: the up gauge matrices, `(n_bands, n_wann, n_kpts)`
- `Udn`: the down gauge matrices, `(n_bands, n_wann, n_kpts)`
"""
function overlap_updn_grad(
        M::AbstractArray{<:Complex, 3},
        Uup::AbstractArray{<:Complex, 3},
        Udn::AbstractArray{<:Complex, 3},
    )
    n_bands, n_wann, n_kpts = size(Uup)

    T = eltype(Uup)
    Mᵂ = zeros(T, n_wann, n_wann)
    GUup = zeros(T, n_bands, n_wann, n_kpts)
    GUdn = zeros(T, n_bands, n_wann, n_kpts)

    for ik in 1:n_kpts
        Uupk = view(Uup, :, :, ik)
        Udnk = view(Udn, :, :, ik)
        Mk = view(M, :, :, ik)
        MUdn = Mk * Udnk
        Mᵂ += Uupk' * MUdn

        GUup[:, :, ik] .= MUdn
        GUdn[:, :, ik] .= Mk' * Uupk
    end

    diagM = diagm(diag(Mᵂ))
    for ik in 1:n_kpts
        view(GUup, :, :, ik) .= view(GUup, :, :, ik) * diagM'
        view(GUdn, :, :, ik) .= view(GUdn, :, :, ik) * diagM
    end
    return GUup ./ n_kpts^2, GUdn ./ n_kpts^2
end

function overlap_updn_grad(model::SpinModel, Uup, Udn)
    return overlap_updn_grad(model.overlaps_updn, Uup, Udn)
end

function omega_updn_grad(model::SpinModel, Uup, Udn)
    # Internal gradient convention is df = 2 Re⟨∇f, dx⟩ (see module docstring
    # on gradient conventions). The minus sign is baked in because
    # `omega_updn` is defined as `n_wann - sum(diag(M))`.
    GUup, GUdn = overlap_updn_grad(model, Uup, Udn)
    return -2 .* GUup, -2 .* GUdn
end

# Canonical-coordinate coupling gradient, taking the ↑↓ overlap directly.
function omega_updn_grad(
        M::AbstractArray{<:Complex, 3},
        Uup::AbstractArray{<:Complex, 3},
        Udn::AbstractArray{<:Complex, 3},
    )
    GUup, GUdn = overlap_updn_grad(M, Uup, Udn)
    return -2 .* GUup, -2 .* GUdn
end
