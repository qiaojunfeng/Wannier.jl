using LinearAlgebra

export spread, center

# -----------------------------------------------------------------------------
# Gradient convention
#
# Internally, every spread-related gradient follows the physics convention
#     df(x) = 2 Re⟨∇f, dx⟩
# i.e. `omega_grad!` and `omega_updn_grad` both return ∇f *in this convention*
# (hence the explicit factor of 2 and the factor of 4 inside `omega_grad!`).
#
# That is a statement about the *derivation*, not about a mismatch to correct:
# the explicit factors of 2 and 4 inside `omega_grad!` are exactly what make the
# emitted gradient the true derivative of the spread in the real-inner-product
# convention Optim.jl / NLSolversBase.jl consume. No rescaling is applied — or
# needed — at the solver boundary, and the finite-difference gradient checks in
# the test suite verify precisely that (analytic == FD, elementwise).
# -----------------------------------------------------------------------------

"""
    imaglog_guided(z, θ)

Branch-stable `Im log z`: evaluates `imaglog(z⋅cis(θ)) − θ`, i.e. picks the
branch of the logarithm closest to `−θ`. With `θ = b ⋅ r_guide` (guiding
centers, as in wannier90) the phases stay continuous across iterations even
when a diagonal overlap approaches the negative real axis. `θ = 0` reproduces
the principal branch.
"""
imaglog_guided(z::Complex, θ::Real) = imaglog(z * cis(θ)) - θ

abstract type AbstractSpread end

@doc raw"""
    struct Spread

The Marzari-Vanderbilt (MV) spread functional.

From MV:
- ``\Omega = \sum_n \langle r^2 \rangle_n - | \langle r \rangle_n |^2``
- ``\langle r \rangle_n = -\frac{1}{N} \sum_{\bm{k},\bm{b}} w_{\bm{b}} \bm{b}
    \Im \log M_{nn}^{\bm{k},\bm{b}}``
- ``\langle r^2 \rangle_n = \frac{1}{N} \sum_{\bm{k},\bm{b}} w_{\bm{b}} \bm{b}
    \left[ \left( 1 - | M_{nn}^{\bm{k},\bm{b}} |^2 \right) +
    \left( \Im \log M_{nn}^{\bm{k},\bm{b}} \right)^2 \right]``

# Fields
- `Ω`: total spread, unit Å²
- `ΩI`: gauge-invarient part, unit Å²
- `ΩOD`: off-diagonal part, unit Å²
- `ΩD`: diagonal part, unit Å²
- `Ωtilde`: Ωtilde = ΩOD + ΩD, unit Å²
- `ω`: Ω of each WF, unit Å², `length(ω) = n_wann`
- `r`: WF center, Cartesian coordinates, unit Å, `3 * n_wann`
"""
struct Spread{T <: Real, C <: Union{Nothing, T}, V <: Union{Nothing, Vector{T}}} <: AbstractSpread
    # Total spread, unit Å², Ω = ΩI + Ωtilde
    Ω::T

    # gauge-invarient part, unit Å²
    ΩI::T

    # off-diagonal part, unit Å²
    ΩOD::T

    # diagonal part, unit Å²
    ΩD::T

    # Ωtilde = ΩOD + ΩD, unit Å². Spelled `tilde` rather than the combining
    # character U+0303 (Ω̃), which is neither typeable nor greppable; the
    # remaining Greek names are the MV paper's own notation (anchor rule).
    Ωtilde::T

    # Ω of each WF, unit Å², length = n_wann
    ω::Vector{T}

    # WF center, Cartesian! coordinates, unit Å, n_wann of Vec3
    r::Vector{Vec3{T}}

    # Optional center-penalty fields. `nothing` when no center penalty applied.
    # additional variables for penalty term
    # Penalty, unit Å²
    Ωc::C

    # Total spread Ωt = Ω + Ωc
    Ωt::C

    # penalty of each WF, unit Å², length = n_wann
    ωc::V

    # total spread of each WF, unit Å², length = n_wann
    # ωt = ω + ωc
    ωt::V
end

function Spread(
        Ω::T, ΩI::T, ΩOD::T, ΩD::T, Ωtilde::T, ω::Vector{T}, r::Vector{Vec3{T}}
    ) where {T <: Real}
    return Spread{T, Nothing, Nothing}(Ω, ΩI, ΩOD, ΩD, Ωtilde, ω, r, nothing, nothing, nothing, nothing)
end

has_center_penalty(Ω::Spread) = Ω.Ωc !== nothing

"""
Preallocated scratch buffers for spread/gradient computation.

`GU` holds the gradient in canonical gauge coordinates, `dΩ/dU*`. It is sized
`(n_bands, n_wannier, n_kpoints)` at construction and never reassigned;
sub-functions that need their own gradient accumulator allocate a separate
buffer and pass it explicitly.
"""
struct Workspace{T}
    X::Array{Complex{T}, 3}
    # Dense gauge-assembly scratch only; `XYLayout` stores just Y's active blocks in x.
    Y::Array{Complex{T}, 3}
    xy::_XYStructure
    U::Array{Complex{T}, 3}
    # gradient in canonical coordinates, dΩ/dU*: n_bands x n_wann x n_kpts
    GU::Array{Complex{T}, 3}
    r::Vector{Vec3{T}}
    # fixed guiding centers for the Im-log branch choice (see
    # `imaglog_guided`); zeros = principal branch (the default). Set them
    # explicitly (e.g. to the trial-orbital centers) when WFs sit far from
    # the origin — like wannier90's `guiding_centres`, they are constants of
    # the run, NOT updated from the current iterate (a self-updating guide
    # would make the objective history-dependent).
    guiding_centers::Vector{Vec3{T}}
    UtMU::Array{Complex{T}, 4}
    MU::Array{Complex{T}, 4}
end

function Workspace(
        bvectors::KspaceStencil{FT}, M::AbstractArray{<:Complex, 4},
        U::AbstractArray{<:Complex, 3},
        frozen::AbstractMatrix{Bool} = falses(size(U, 1), size(U, 3)),
    ) where {FT}
    n_kpts = size(M, 4)
    n_bands = size(M, 1)
    n_wann = size(U, 2)
    n_bvecs = size(M, 3)

    X = zeros(Complex{FT}, n_wann, n_wann, n_kpts)
    Y = zeros(Complex{FT}, n_bands, n_wann, n_kpts)
    xy = _xy_structure(frozen, n_wann)
    _initialize_compact_y!(Y, xy)
    Ucopy = zeros(Complex{FT}, n_bands, n_wann, n_kpts)
    GU = zeros(Complex{FT}, n_bands, n_wann, n_kpts)
    r = zeros(Vec3{FT}, n_wann)

    MU = zeros(Complex{FT}, n_bands, n_wann, n_bvecs, n_kpts)
    UtMU = zeros(Complex{FT}, n_wann, n_wann, n_bvecs, n_kpts)

    return Workspace(X, Y, xy, Ucopy, GU, r, zeros(Vec3{FT}, n_wann), UtMU, MU)
end

Workspace(model::Model) =
    Workspace(model.kstencil, model.overlaps, model.gauges, model.frozen_bands)

n_bands(w::Workspace) = size(w.GU, 1)
n_wannier(w::Workspace) = size(w.GU, 2)
n_kpoints(w::Workspace) = size(w.GU, 3)

function _alloc_mu_utmu_packed(::Type{FT}, n_kpts, n_bvecs, n_bands, n_wann) where {FT}
    MU = zeros(Complex{FT}, n_bands, n_wann, n_bvecs, n_kpts)
    UtMU = zeros(Complex{FT}, n_wann, n_wann, n_bvecs, n_kpts)
    return MU, UtMU
end

center_penalty(r0, λ) = (r, n) -> (r - λ * (r - r0[n]))

"""
    omega_center(bvectors, M, U, r0, λ)

Compute WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
- `r0`: `3 * n_wann`, WF centers in cartesian coordinates
- `λ`: penalty strength
"""
function omega_center(args...; kwargs...)
    Ω = spread(args...)
    return omega_center(Ω; kwargs...)
end

function omega_center(Ω::Spread; r0::Vector{Vec3{T}}, λ::T) where {T <: Real}
    ωc = λ .* map(i -> (t = Ω.r[i] - r0[i]; sum(t .^ 2)), 1:length(r0))
    ωt = Ω.ω + ωc
    Ωc = sum(ωc)
    Ωt = Ω.Ω + Ωc
    return Spread(Ω.Ω, Ω.ΩI, Ω.ΩOD, Ω.ΩD, Ω.Ωtilde, Ω.ω, Ω.r, Ωc, Ωt, ωc, ωt)
end

function omega!(
        r::Vector{<:Vec3{FT}},
        UtMU::AbstractArray{<:Complex, 4},
        MU::AbstractArray{<:Complex, 4},
        bvectors::KspaceStencil{FT},
        M;
        guide::AbstractVector{<:Vec3} = zeros(Vec3{FT}, length(r)),
    ) where {FT <: Real}
    rg = guide
    fill!(r, zero(eltype(r)))

    nw = length(r)
    nk = size(UtMU, 4)
    n_bvecs = size(UtMU, 3)

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = reciprocal_lattice(bvectors)
    kpoints = bvectors.kpoints

    r² = zeros(FT, nw)

    ΩI::FT = 0.0
    ΩOD::FT = 0.0
    ΩD::FT = 0.0

    @inbounds for ik in 1:nk
        for ib in 1:n_bvecs
            ikpb = kpb_k[ib, ik]
            Nkb = view(UtMU, :, :, ib, ik)
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ib, ik] - kpoints[ik])

            wᵇ = wb[ib]
            wb_b = wᵇ * b

            ts = zero(FT)
            ts2 = zero(FT)
            for i in axes(Nkb, 2)
                for j in axes(Nkb, 1)
                    nt = Nkb[j, i]
                    a2 = abs2(nt)
                    ts += a2

                    if i == j
                        imlogN = imaglog_guided(nt, b ⋅ rg[i])

                        r[i] -= imlogN * wb_b
                        r²[i] += wᵇ * (1 - a2 + imlogN^2)
                    else
                        ts2 += a2
                    end
                end
            end

            ΩI += wᵇ * (nw - ts)
            ΩOD += wᵇ * ts2
        end
    end

    # in place: the caller's workspace keeps the true centers (they seed the
    # Im-log branch guides of the next evaluation)
    r .= r ./ nk
    r² /= nk
    ΩI /= nk
    ΩOD /= nk

    @inbounds for ik in 1:nk
        for ib in 1:n_bvecs
            ikpb = kpb_k[ib, ik]
            Nᵏᵇ = view(UtMU, :, :, ib, ik)
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ib, ik] - kpoints[ik])
            wᵇ = wb[ib]

            for n in 1:nw
                ΩD += wᵇ * (-imaglog_guided(Nᵏᵇ[n, n], b ⋅ rg[n]) - b' * r[n])^2
            end
        end
    end

    ΩD /= nk
    Ωtilde = ΩOD + ΩD
    ω = r² - map(x -> sum(abs.(x .^ 2)), r)
    Ω = ΩI + Ωtilde

    return Spread(Ω, ΩI, ΩOD, ΩD, Ωtilde, ω, r)
end

function omega!(cache::Workspace, bvectors::KspaceStencil{FT}, M) where {FT <: Real}
    return omega!(cache.r, cache.UtMU, cache.MU, bvectors, M; guide = cache.guiding_centers)
end

"""
    spread(model, [U])
    spread(bvectors, M, U)

Compute WF spread for a [`Model`](@ref), potentially for a given gauge `U`, or by explicitely giving
`bvectors` and `M`.
In case of the first `bvectors = model.bvectors` and `M = model.overlaps_updn`.
"""
spread(model::Model) = spread(model, model.gauges)
spread(model::Model, gauges) = spread(model.kstencil, model.overlaps, gauges)

function spread(bvectors::KspaceStencil, M::AbstractArray{<:Complex, 4}, U::AbstractArray{<:Complex, 3})
    n_kpts = size(M, 4)
    n_bvecs = size(M, 3)
    n_bands = size(U, 1)
    n_wann = size(U, 2)
    FT = real(eltype(M))
    MU, UtMU = _alloc_mu_utmu_packed(FT, n_kpts, n_bvecs, n_bands, n_wann)
    r = zeros(Vec3{FT}, n_wann)

    compute_MU_UtMU!(MU, UtMU, bvectors, M, U)
    return omega!(r, UtMU, MU, bvectors, M)
end

function Base.show(io::IO, ::MIME"text/plain", Ω::Spread)
    n_wann = length(Ω.ω)
    if has_center_penalty(Ω)
        println(io, "  WF     center [rx, ry, rz]/Å              spread/Å²  ω  ωc  ωt")
        for i in 1:n_wann
            @printf(
                io,
                "%4d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n",
                i, Ω.r[i]..., Ω.ω[i], Ω.ωc[i], Ω.ωt[i],
            )
        end
        @printf(io, "Sum spread: Ωt = Ω + Ωc, Ω = ΩI + Ωtilde, Ωtilde = ΩOD + ΩD\n")
        @printf(io, "   Ωt  = %11.5f\n", Ω.Ωt)
        @printf(io, "   Ωc  = %11.5f\n", Ω.Ωc)
        @printf(io, "   Ω   = %11.5f\n", Ω.Ω)
        @printf(io, "   ΩI  = %11.5f\n", Ω.ΩI)
        @printf(io, "   ΩOD = %11.5f\n", Ω.ΩOD)
        @printf(io, "   ΩD  = %11.5f\n", Ω.ΩD)
        return @printf(io, "   Ωtilde   = %11.5f", Ω.Ωtilde)
    else
        println(io, "  WF     center [rx, ry, rz]/Å              spread/Å²")
        for i in 1:n_wann
            @printf(io, "%4d %11.5f %11.5f %11.5f %11.5f\n", i, Ω.r[i]..., Ω.ω[i])
        end
        @printf(io, "Sum spread: Ω = ΩI + Ωtilde, Ωtilde = ΩOD + ΩD\n")
        @printf(io, "   ΩI  = %11.5f\n", Ω.ΩI)
        @printf(io, "   Ωtilde   = %11.5f\n", Ω.Ωtilde)
        @printf(io, "   ΩOD = %11.5f\n", Ω.ΩOD)
        @printf(io, "   ΩD  = %11.5f\n", Ω.ΩD)
        return @printf(io, "   Ω   = %11.5f\n", Ω.Ω)
    end
end

omega_grad!(cache::Workspace, bvectors, M) = omega_grad!((r, _) -> r, cache, bvectors, M)

"""Gradient with an externally-provided buffer `G`; leaves `cache.GU` untouched."""
omega_grad!(G::AbstractArray{<:Complex, 3}, cache::Workspace, bvectors, M) =
    omega_grad!((r, _) -> r, G, cache.r, cache.UtMU, cache.MU, bvectors, M; rg = cache.guiding_centers)

omega_grad!(penalty::Function, G::AbstractArray{<:Complex, 3}, cache::Workspace, bvectors, M) =
    omega_grad!(penalty, G, cache.r, cache.UtMU, cache.MU, bvectors, M; rg = cache.guiding_centers)

function omega_grad!(
        penalty::Function,
        G::AbstractArray{<:Complex, 3},
        r,
        UtMU::AbstractArray{<:Complex, 4},
        MU::AbstractArray{<:Complex, 4},
        bvectors,
        M;
        rg::AbstractVector{<:Vec3} = zeros(eltype(r), length(r)),
    )
    fill!(G, 0)

    n_bands, n_wann, n_kpts = size(G)
    scratch = zeros(eltype(G), n_bands, n_wann)

    n_bvecs = size(UtMU, 3)

    center!(r, UtMU, bvectors; guide = rg)

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    @inbounds for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ib, ik]
            MUkb = view(MU, :, :, ib, ik)
            Nkb = view(UtMU, :, :, ib, ik)
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ib, ik] - kpoints[ik])

            wᵇ = wb[ib]

            for n in 1:n_wann
                nn = Nkb[n, n]
                q = imaglog_guided(nn, b ⋅ rg[n]) + penalty(r[n], n) ⋅ b
                t = -im * q / nn
                cnn = conj(nn)
                for m in 1:n_bands
                    scratch[m, n] = MUkb[m, n] * (t - cnn)
                end
            end

            view(G, :, :, ik) .+= 4 .* wᵇ .* scratch
        end
    end

    G ./= n_kpts

    return G
end

function omega_grad!(penalty::Function, cache::Workspace{T}, bvectors, M) where {T}
    return omega_grad!(penalty, cache.GU, cache.r, cache.UtMU, cache.MU, bvectors, M; rg = cache.guiding_centers)
end

"""
    omega_grad(bvectors, M, U, r)

Compute gradient of WF spread.

Size of output `dΩ/dU` = `n_bands * n_wann * n_kpts`.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
- `r`: `3 * n_wann`, the current WF centers in cartesian coordinates
"""
function omega_grad(penalty::Function, bvectors::KspaceStencil, M::AbstractArray{<:Complex, 4}, U::AbstractArray{<:Complex, 3})
    n_kpts = size(M, 4)
    n_bvecs = size(M, 3)
    n_bands = size(U, 1)
    n_wann = size(U, 2)
    FT = real(eltype(M))
    MU, UtMU = _alloc_mu_utmu_packed(FT, n_kpts, n_bvecs, n_bands, n_wann)
    r = zeros(Vec3{FT}, n_wann)
    G = zeros(Complex{FT}, n_bands, n_wann, n_kpts)

    compute_MU_UtMU!(MU, UtMU, bvectors, M, U)
    return omega_grad!(penalty, G, r, UtMU, MU, bvectors, M)
end
omega_grad(bvectors::KspaceStencil, M, U) = omega_grad((r, _) -> r, bvectors, M, U)

"""
    center(bvectors, M, U)

Compute WF center in reciprocal space.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
function center(bvectors::KspaceStencil, M::AbstractArray{<:Complex, 4}, U::AbstractArray{<:Complex, 3})
    n_kpts = size(M, 4)
    n_bvecs = size(M, 3)
    n_bands = size(U, 1)
    n_wann = size(U, 2)
    FT = real(eltype(M))
    MU, UtMU = _alloc_mu_utmu_packed(FT, n_kpts, n_bvecs, n_bands, n_wann)
    r = zeros(Vec3{FT}, n_wann)

    compute_MU_UtMU!(MU, UtMU, bvectors, M, U)
    return center!(r, UtMU, bvectors)
end

function center!(
        r::Vector{<:Vec3}, UtMU::AbstractArray{<:Complex, 4}, bvectors;
        guide::AbstractVector{<:Vec3} = zeros(eltype(r), length(r)),
    )
    rg = guide
    fill!(r, zero(eltype(r)))
    n_wann = length(r)

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = reciprocal_lattice(bvectors)
    kpoints = bvectors.kpoints

    nk = size(UtMU, 4)
    nb = size(UtMU, 3)

    @inbounds for ik in 1:nk
        k = kpoints[ik]
        for ib in 1:nb
            Nb = view(UtMU, :, :, ib, ik)
            ikpb = kpb_k[ib, ik]

            b = recip_lattice * (kpoints[ikpb] + kpb_G[ib, ik] - k)
            w = wb[ib]

            for n in 1:n_wann
                fac = w * imaglog_guided(Nb[n, n], b ⋅ rg[n])
                r[n] -= b * fac
            end
        end
    end

    r ./= nk

    return r
end

"""
    center(model)

Compute WF center in reciprocal space for `Model`.
"""
center(model::Model) = center(model.kstencil, model.overlaps, model.gauges)

"""
    center(model, U)

Compute WF center in reciprocal space for `Model` with given `U` gauge.

# Arguments
- `model`: the `Model`
- `U`: `n_bands × n_wann × n_kpts` array
"""
function center(model::Model, U::AbstractArray{<:Complex, 3})
    return center(model.kstencil, model.overlaps, U)
end

"""
    position_operator(bvectors, M, U)

Compute WF postion operator matrix in reciprocal space.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
@views function position_operator(
        bvectors::KspaceStencil{FT},
        M::AbstractArray{Complex{FT}, 4},
        U::AbstractArray{Complex{FT}, 3},
    ) where {FT <: Real}
    n_bands = size(U, 1)
    n_wann = size(U, 2)
    n_kpts = size(U, 3)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # along x, y, z directions
    R = zeros(Vec3{Complex{FT}}, n_wann, n_wann)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ib, ik]

            Nᵏᵇ .= view(U, :, :, ik)' * M[:, :, ib, ik] * view(U, :, :, ikpb)
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ib, ik] - kpoints[ik])

            wᵇ = wb[ib]

            for m in 1:n_wann
                for n in 1:n_wann
                    R[m, n] += wᵇ * Nᵏᵇ[m, n] * b

                    if m == n
                        R[m, n] -= wᵇ * b
                    end
                end
            end
        end
    end

    R /= -im * n_kpts

    return R
end

"""
    position_operator(model)

Compute WF postion operator matrix in reciprocal space for `Model`.
"""
position_operator(model::Model) = position_operator(model.kstencil, model.overlaps, model.gauges)

"""
    position_operator(model, U)

Compute WF postion operator matrix in reciprocal space for `Model` with given `U` gauge.

# Arguments
- `model`: the `Model`
- `U`: `n_bands × n_wann × n_kpts` array
"""
function position_operator(model::Model, U::AbstractArray{<:Complex, 3})
    return position_operator(model.kstencil, model.overlaps, U)
end

"""
    $(SIGNATURES)

Wannier-gauge Berry connection in kspace, WYSV Eq. 44 or MV Eq. C14

# Keyword arguments
- `imlog_diag`: use imaginary part of logrithm for diagonal elements,
    MV1997 Eq. 31. wannier90 default is true.
"""
function compute_berry_connection_kspace(
        kstencil::KspaceStencil,
        overlaps::AbstractArray{<:Complex, 4},
        gauges::AbstractArray{<:Complex, 3};
        imlog_diag::Bool = true,
        force_hermiticity::Bool = default_w90_berry_position_force_hermiticity(),
    )
    nkpts = size(gauges, 3)
    @assert nkpts > 0 "empty gauges"
    nwann = size(gauges, 2)
    nbvecs = size(overlaps, 3)

    wb = kstencil.bweights
    recip_lattice = reciprocal_lattice(kstencil)
    kpoints = kstencil.kpoints

    # Aᵂ can be indexed by Aᵂ[ik][m, n][α], where
    # - ik: kpoint index
    # - m, n: WF indices
    # - α ∈ {1, 2, 3} for x, y, z directions
    Aᵂ = map(1:nkpts) do ik
        Uₖ = view(gauges, :, :, ik)
        Aᵂₖ = zeros(Vec3{eltype(Uₖ)}, nwann, nwann)
        for ib in 1:nbvecs
            ik2 = kstencil.kpb_k[ib, ik]
            Mᴴ = view(overlaps, :, :, ib, ik)
            Uₖ₂ = view(gauges, :, :, ik2)
            Mᵂ = Uₖ' * Mᴴ * Uₖ₂
            G = kstencil.kpb_G[ib, ik]
            # b isa Vec3, along x, y, z directions
            b = recip_lattice * (kpoints[ik2] + G - kpoints[ik])
            if imlog_diag
                N = im * (Mᵂ - I)
                N[diagind(N)] .= -imaglog.(diag(Mᵂ))
                Aᵂₖ .+= wb[ib] .* Ref(b) .* N
            else
                Aᵂₖ .+= im * wb[ib] .* Ref(b) .* (Mᵂ - I)
            end
        end
        if force_hermiticity
            # cannot use adjoint which recursively transpose the inner Vec3
            Aᵂₖ = (Aᵂₖ + conj(permutedims(Aᵂₖ))) / 2
        end
        return Aᵂₖ
    end
    return Aᵂ
end

@inline function compute_berry_connection_kspace(
        model::Model, gauges::AbstractArray{<:Complex, 3} = model.gauges; kwargs...
    )
    return compute_berry_connection_kspace(
        model.kstencil, model.overlaps, gauges; kwargs...
    )
end

"""
    omega_local(bvectors, M, U)

Local part of the contribution to `r^2`.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * n_bvecs * n_kpts` overlap array
- `U`: `n_bands × n_wann × n_kpts` array
"""
function omega_local(
        bvectors::KspaceStencil{FT},
        M::AbstractArray{Complex{FT}, 4},
        U::AbstractArray{Complex{FT}, 3},
    ) where {FT <: Real}
    n_bands = size(U, 1)
    n_wann = size(U, 2)
    n_kpts = size(U, 3)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    wb = bvectors.bweights

    loc = zeros(FT, n_kpts)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ib, ik]
            Nᵏᵇ .= view(U, :, :, ik)' * M[:, :, ib, ik] * view(U, :, :, ikpb)

            for n in 1:n_wann
                loc[ik] += wb[ib] * (1 - abs(Nᵏᵇ[n, n])^2 + imaglog(Nᵏᵇ[n, n])^2)
            end
        end
    end

    return loc
end
