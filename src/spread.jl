using LinearAlgebra

export omega, omega_grad, center, SpreadPenalty, CenterSpreadPenalty

abstract type AbstractSpread end
abstract type AbstractPenalty end

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
- `Ω̃`: Ω̃ = ΩOD + ΩD, unit Å²
- `ω`: Ω of each WF, unit Å², `length(ω) = n_wann`
- `r`: WF center, Cartesian coordinates, unit Å, `3 * n_wann`
"""
struct Spread{T<:Real} <: AbstractSpread
    # Total spread, unit Å², Ω = ΩI + Ω̃
    Ω::T

    # gauge-invarient part, unit Å²
    ΩI::T

    # off-diagonal part, unit Å²
    ΩOD::T

    # diagonal part, unit Å²
    ΩD::T

    # Ω̃ = ΩOD + ΩD, unit Å²
    Ω̃::T

    # Ω of each WF, unit Å², length = n_wann
    ω::Vector{T}

    # WF center, Cartesian! coordinates, unit Å, n_wann of Vec3
    r::Vector{Vec3{T}}

    # frozen_weight::T
    # fix_centers :: Array{Float64,2} #3 x nwannier
end

# TODO refactor, this is a copy-paste of `Spread` :-(
"""
    struct SpreadCenter

A `struct` containing both `Spread` and WF center penalty.
"""
struct SpreadCenter{T} <: AbstractSpread
    # Total spread, unit Å², Ω = ΩI + Ω̃
    Ω::T

    # gauge-invarient part, unit Å²
    ΩI::T

    # off-diagonal part, unit Å²
    ΩOD::T

    # diagonal part, unit Å²
    ΩD::T

    # Ω̃ = ΩOD + ΩD, unit Å²
    Ω̃::T

    # Ω of each WF, unit Å², length = n_wann
    ω::Vector{T}

    # WF center, Cartesian! coordinates, unit Å, 3 * n_wann
    r::Vector{Vec3{T}}

    # additional variables for penalty term
    # Penalty, unit Å²
    Ωc::T

    # Total spread Ωt = Ω + Ωc
    Ωt::T

    # penalty of each WF, unit Å², length = n_wann
    ωc::Vector{T}

    # total spread of each WF, unit Å², length = n_wann
    # ωt = ω + ωc
    ωt::Vector{T}
end

function Base.show(io::IO, ::MIME"text/plain", Ω::SpreadCenter)
    println(io, "  WF     center [rx, ry, rz]/Å              spread/Å²  ω  ωc  ωt")

    n_wann = length(Ω.ω)

    for i in 1:n_wann
        @printf(
            io,
            "%4d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n",
            i,
            Ω.r[i]...,
            Ω.ω[i],
            Ω.ωc[i],
            Ω.ωt[i]
        )
    end

    @printf(io, "Sum spread: Ωt = Ω + Ωc, Ω = ΩI + Ω̃, Ω̃ = ΩOD + ΩD\n")
    @printf(io, "   Ωt  = %11.5f\n", Ω.Ωt)
    @printf(io, "   Ωc  = %11.5f\n", Ω.Ωc)
    @printf(io, "   Ω   = %11.5f\n", Ω.Ω)
    @printf(io, "   ΩI  = %11.5f\n", Ω.ΩI)
    @printf(io, "   ΩOD = %11.5f\n", Ω.ΩOD)
    @printf(io, "   ΩD  = %11.5f\n", Ω.ΩD)
    @printf(io, "   Ω̃   = %11.5f", Ω.Ω̃)
end

#TODO a bit generic
mutable struct Cache{T}
    X::Vector{Matrix{Complex{T}}}
    Y::Vector{Matrix{Complex{T}}}
    U::Vector{Matrix{Complex{T}}}
    # n_bands x n_wann x n_kpts
    G::Array{Complex{T},3}
    r::Vector{Vec3{T}}
    UtMU::Vector{Vector{Matrix{Complex{T}}}}
    MU::Vector{Vector{Matrix{Complex{T}}}}
end

# TODO remove 1st arg?
function Cache(
    bvectors::KspaceStencil{FT}, M::Vector{Vector{Matrix{Complex{FT}}}}, U
) where {FT}
    n_kpts = length(M)
    n_bands, n_wann = U isa Vector ? size(U[1]) : size.((U,), (1, 2))
    n_bvecs = length(M[1])

    X = [zeros(Complex{FT}, n_wann, n_wann) for i in 1:n_kpts]
    Y = [zeros(Complex{FT}, n_bands, n_wann) for i in 1:n_kpts]
    U = [zeros(Complex{FT}, n_bands, n_wann) for i in 1:n_kpts]
    G = zeros(Complex{FT}, n_bands, n_wann, n_kpts)
    r = zeros(Vec3{FT}, n_wann)

    UtMU = [[zeros(Complex{FT}, n_wann, n_wann) for ib in 1:n_bvecs] for i in 1:n_kpts]
    MU = [[zeros(Complex{FT}, n_bands, n_wann) for ib in 1:n_bvecs] for i in 1:n_kpts]

    return Cache(X, Y, U, G, r, UtMU, MU)
end

Cache(model::Model) = Cache(model.kstencil, model.overlaps, model.gauges)

n_bands(c::Cache) = size(c.G, 1)
n_wann(c::Cache) = size(c.G, 2)
n_kpts(c::Cache) = size(c.G, 3)

function compute_MU_UtMU!(MU, UtMU, bvectors::KspaceStencil, M, U::Vector)
    kpb_k = bvectors.kpb_k
    n_bvecs = length(kpb_k[1])

    @inbounds for ik in 1:length(U)
        Ut = U[ik]'
        for ib in 1:n_bvecs
            MUkb = MU[ik][ib]
            ikpb = kpb_k[ik][ib]
            mul!(MUkb, M[ik][ib], U[ikpb])
            mul!(UtMU[ik][ib], Ut, MUkb)
        end
    end
    return MU, UtMU
end

function compute_MU_UtMU!(MU, UtMU, bvectors::KspaceStencil, M, U::Array)
    kpb_k = bvectors.kpb_k
    n_bvecs = length(kpb_k[1])

    @inbounds for ik in axes(U, 3)
        Ut = view(U, :, :, ik)'
        for ib in 1:n_bvecs
            MUkb = MU[ik][ib]
            ikpb = kpb_k[ik][ib]
            Ukpb = view(U, :, :, ikpb)
            mul!(MUkb, M[ik][ib], Ukpb)
            mul!(UtMU[ik][ib], Ut, MUkb)
        end
    end
    return MU, UtMU
end
function compute_MU_UtMU!(cache::Cache, bvectors::KspaceStencil, M, U)
    return compute_MU_UtMU!(cache.MU, cache.UtMU, bvectors, M, U)
end

"""
Standard penalty for minimizing the total spread.
"""
struct SpreadPenalty <: AbstractPenalty end

omega!(::SpreadPenalty, args...) = omega!(args...)
omega(::SpreadPenalty, args...) = omega(args...)

omega_grad!(::SpreadPenalty, args...) = omega_grad!(args...)
omega_grad(::SpreadPenalty, args...) = omega_grad(args...)

"""
Penalty for minimizing the spread as well as maximizing the "closeness" to the atoms.
"""
struct CenterSpreadPenalty{T} <: AbstractPenalty
    r₀::Vector{Vec3{T}}
    λ::T
end

# TODO probably omega should always just return a value not the spread struct
omega!(p::CenterSpreadPenalty, args...) = (Ω=omega_center(omega!(args...); p.r₀, p.λ).Ωt,)
omega(p::CenterSpreadPenalty, args...) = omega_center(omega(args...); p.r₀, p.λ)

function omega_grad!(p::CenterSpreadPenalty, args...)
    return omega_grad!(center_penalty(p.r₀, p.λ), args...)
end

center_penalty(r₀, λ) = (r, n) -> (r - λ * (r - r₀[n]))

"""
    omega_center(bvectors, M, U, r₀, λ)

Compute WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
- `r₀`: `3 * n_wann`, WF centers in cartesian coordinates
- `λ`: penalty strength
"""
function omega_center(args...; kwargs...)
    Ω = omega(args...)
    return omega_center(Ω; kwargs...)
end

function omega_center(Ω::Spread; r₀::Vector{Vec3{T}}, λ::T) where {T<:Real}
    ωc = λ .* map(i -> (t = Ω.r[i] - r₀[i]; sum(t .^ 2)), 1:length(r₀))
    ωt = Ω.ω + ωc
    Ωc = sum(ωc)
    Ωt = Ω.Ω + Ωc
    return SpreadCenter(Ω.Ω, Ω.ΩI, Ω.ΩOD, Ω.ΩD, Ω.Ω̃, Ω.ω, Ω.r, Ωc, Ωt, ωc, ωt)
end

function omega!(cache::Cache, bvectors::KspaceStencil{FT}, M) where {FT<:Real}
    r = cache.r
    fill!(r, zero(eltype(r)))

    UtMU = cache.UtMU
    MU = cache.MU

    nw = n_wann(cache)
    nk = n_kpts(cache)

    n_bvecs = length(M[1])

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = reciprocal_lattice(bvectors)
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    r² = zeros(FT, nw)

    ΩI::FT = 0.0
    ΩOD::FT = 0.0
    ΩD::FT = 0.0

    for ik in 1:nk
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])
        MUk = MU[ik]
        Nk = UtMU[ik]

        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            MUkb = MUk[ib]
            Nkb = Nk[ib]
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ik][ib] - kpoints[ik])

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
                        imlogN = imaglog(nt)

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

    r = map(x -> x ./ nk, r)
    r² /= nk
    ΩI /= nk
    ΩOD /= nk
    # w_froz /= n_kpts

    # ΩD requires r, so we need different loops
    # However, since ΩD = Ω - ΩI - ΩOD, we can skip these loops
    # for ik in 1:n_kpts
    #     for ib in 1:n_bvecs
    #         ikpb = kpb_k[ib, ik]
    #         Nᵏᵇ .= U[:, :, ik]' * M[:, :, ib, ik] * U[:, :, ikpb]
    #         b .= recip_lattice * (kpoints[:, ikpb] + kpb_G[:, ib, ik] - kpoints[:, ik])
    #         wᵇ = wb[ib]

    #         for n in 1:n_wann
    #             ΩD += wᵇ * (-imaglog(Nᵏᵇ[n, n]) - b' * r[:, n])^2
    #         end
    #     end
    # end
    # ΩD /= n_kpts
    # Ω̃ = ΩOD + ΩD

    # @debug "Spread" r r²'
    # @debug "Spread" ΩI ΩOD ΩD

    # Ω of each WF
    ω = r² - map(x -> sum(abs.(x .^ 2)), r)
    # total Ω
    Ω = sum(ω)
    # Ω += w_froz
    Ω̃ = Ω - ΩI
    ΩD = Ω̃ - ΩOD

    return Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r)
    # return Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r, w_froz)
end

"""
    omega(model, [U])
    omega(bvectors, M, U)

Compute WF spread for a [`Model`](@ref), potentially for a given gauge `U`, or by explicitely giving
`bvectors` and `M`.
In case of the first `bvectors = model.bvectors` and `M = model.M`.
"""
omega(model::Model) = omega(model, model.gauges)
omega(model::Model, gauges) = omega(model.kstencil, model.overlaps, gauges)
function omega(bvectors::KspaceStencil, M, X, Y)
    U = X_Y_to_U(X, Y)
    return omega(bvectors, M, U)
end

function omega(bvectors::KspaceStencil, M, U)
    cache = Cache(bvectors, M, U)
    compute_MU_UtMU!(cache, bvectors, M, U)
    return omega!(cache, bvectors, M)
end

function Base.show(io::IO, ::MIME"text/plain", Ω::Spread)
    println(io, "  WF     center [rx, ry, rz]/Å              spread/Å²")

    n_wann = length(Ω.ω)
    for i in 1:n_wann
        @printf(io, "%4d %11.5f %11.5f %11.5f %11.5f\n", i, Ω.r[i]..., Ω.ω[i])
    end

    @printf(io, "Sum spread: Ω = ΩI + Ω̃, Ω̃ = ΩOD + ΩD\n")
    @printf(io, "   ΩI  = %11.5f\n", Ω.ΩI)
    @printf(io, "   Ω̃   = %11.5f\n", Ω.Ω̃)
    @printf(io, "   ΩOD = %11.5f\n", Ω.ΩOD)
    @printf(io, "   ΩD  = %11.5f\n", Ω.ΩD)
    @printf(io, "   Ω   = %11.5f\n", Ω.Ω)
end
omega_grad!(cache::Cache, bvectors, M) = omega_grad!((r, _) -> r, cache, bvectors, M)
function omega_grad!(penalty::Function, cache::Cache{T}, bvectors, M) where {T}
    # This mutates cache.G and cache.Mkb
    G = cache.G
    fill!(G, 0)
    r = cache.r
    UtMU = cache.UtMU
    MU = cache.MU

    n_bands, n_wann, n_kpts = size(G)

    n_bvecs = length(M[1])

    center!(r, UtMU, bvectors)

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    @inbounds for ik in 1:n_kpts
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])
        # G[1:n_froz, :, ik] = -2 * μ * U[1:n_froz, :, ik]
        MUk = MU[ik]
        Nk = UtMU[ik]
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            MUkb = MUk[ib]
            Nkb = Nk[ib]
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ik][ib] - kpoints[ik])

            wᵇ = wb[ib]

            # MV way
            # fA(B) = (B - B') / 2
            # fS(B) = (B + B') / (2 * im)
            # q = imaglog.(diag(Nᵏᵇ)) + r' * b
            # for m = 1:n_wann, n = 1:n_wann
            #     R[m, n] = Nᵏᵇ[m, n] * conj(Nᵏᵇ[n, n])
            #     T[m, n] = Nᵏᵇ[m, n] / Nᵏᵇ[n, n] * q[n]
            # end
            # G[:, :, ik] += 4 * wᵇ * (fA(R) .- fS(T))

            for n in 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                nn = Nkb[n, n]

                # TODO: This check can be done somewherhe else adds 12% of time
                # if abs(nn) < 1e-10
                #     display(Nk[ib])
                #     println()
                #     error("Nᵏᵇ too small! $ik -> $ikpb")
                # end

                q = imaglog(nn) + penalty(r[n], n) ⋅ b

                t = -im * q / nn

                cnn = conj(nn)
                for m in 1:n_bands
                    # T[m, n] = -im * MUᵏᵇ[m, n] / (Nᵏᵇ[n, n]) * q[n]
                    # TODO: this changes MU[ik][ib], so the cache.MU is not
                    # the mathematically correct MU anymore.
                    # should refactor this to make sure quantities in cache
                    # is mathematically correct
                    MUkb[m, n] *= (t - cnn)
                end
            end

            view(G, :, :, ik) .+= 4 .* wᵇ .* MUk[ib]
        end
    end

    G ./= n_kpts

    return G
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
function omega_grad(penalty::Function, bvectors::KspaceStencil, M, U)
    cache = Cache(bvectors, M, U)
    compute_MU_UtMU!(cache, bvectors, M, U)
    return omega_grad!(penalty, cache, bvectors, M)
end
omega_grad(bvectors::KspaceStencil, M, U) = omega_grad((r, _) -> r, bvectors, M, U)

function omega_grad(penalty::Function, bvectors::KspaceStencil, M, X, Y, frozen)
    U = X_Y_to_U(X, Y)
    G = omega_grad(penalty, bvectors, M, U)
    return GU_to_GX_GY(G, X, Y, frozen)
end
function omega_grad(bvectors::KspaceStencil, M, X, Y, frozen)
    return omega_grad((r, _) -> r, bvectors, M, X, Y, frozen)
end

"""
    omega_local(bvectors, M, U)

Local part of the contribution to `r^2`.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
function omega_local(
    bvectors::KspaceStencil{FT}, M::Vector, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = length(M[1])

    kpb_k = bvectors.kpb_k
    wb = bvectors.bweights

    loc = zeros(FT, n_kpts)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            Nᵏᵇ .= U[ik]' * M[ik][ib] * U[ikpb]

            for n in 1:n_wann
                loc[ik] += wb[ib] * (1 - abs(Nᵏᵇ[n, n])^2 + imaglog(Nᵏᵇ[n, n])^2)
            end
        end
    end

    return loc
end

"""
    center(bvectors, M, U)

Compute WF center in reciprocal space.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
function center(bvectors::KspaceStencil, M, U)
    cache = Cache(bvectors, M, U)
    compute_MU_UtMU!(cache, bvectors, M, U)
    return center!(cache.r, cache.UtMU, bvectors)
end
function center!(r::Vector{<:Vec3}, UtMU, bvectors)
    fill!(r, zero(eltype(r)))
    n_wann = length(r)

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = reciprocal_lattice(bvectors)
    kpoints = bvectors.kpoints

    @inbounds for (ik, Nk) in enumerate(UtMU)
        k = kpoints[ik]
        for (ib, Nb) in enumerate(Nk)
            ikpb = kpb_k[ik][ib]

            b = recip_lattice * (kpoints[ikpb] + kpb_G[ik][ib] - k)

            w = wb[ib]

            for n in 1:n_wann
                fac = w * imaglog(Nb[n, n])
                r[n] -= b * fac
            end
        end
    end

    r ./= length(UtMU)

    return r
end

"""
    center(bvectors, M, U)

Compute WF center in reciprocal space.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
function center(
    bvectors::KspaceStencil{FT}, M::Vector, U::Array{Complex{FT},3}
) where {FT<:Real}
    n_bands, n_wann, n_kpts = size(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    r = zeros(Vec3{FT}, n_wann)
    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    cache = zeros(Complex{FT}, n_bands, n_wann)
    # M_ = map(ik -> map(ib -> , 1:n_bvecs), 1:n_kpts)

    @inbounds @views for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            mul!(cache, M[ik][ib], U[:, :, ikpb])
            mul!(Nᵏᵇ, U[:, :, ik]', cache)
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ik][ib] - kpoints[ik])

            w = wb[ib]

            for n in 1:n_wann
                fac = w * imaglog(Nᵏᵇ[n, n])
                r[n] -= b * fac
            end
        end
    end

    r ./= n_kpts

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
- `U`: `n_wann * n_wann * n_kpts` array
"""
function center(model::Model, U::Vector{Matrix{T}}) where {T<:Number}
    return center(model.kstencil, model.overlaps, U)
end

"""
    position_op(bvectors, M, U)

Compute WF postion operator matrix in reciprocal space.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
@views function position_op(
    bvectors::KspaceStencil{FT}, M::Vector, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = length(M[1])

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
            ikpb = kpb_k[ik][ib]

            Nᵏᵇ .= U[ik]' * M[ik][ib] * U[ikpb]
            b = recip_lattice * (kpoints[ikpb] + kpb_G[ik][ib] - kpoints[ik])

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
    position_op(model)

Compute WF postion operator matrix in reciprocal space for `Model`.
"""
position_op(model::Model) = position_op(model.kstencil, model.overlaps, model.gauges)

"""
    position_op(model, U)

Compute WF postion operator matrix in reciprocal space for `Model` with given `U` gauge.

# Arguments
- `model`: the `Model`
- `U`: `n_wann * n_wann * n_kpts` array
"""
function position_op(model::Model, U::Vector{Matrix{T}}) where {T<:Number}
    return position_op(model.kstencil, model.overlaps, U)
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
    overlaps::AbstractVector,
    gauges::AbstractVector;
    imlog_diag::Bool=true,
    force_hermiticity::Bool=default_w90_berry_position_force_hermiticity(),
)
    nkpts = length(gauges)
    @assert nkpts > 0 "empty gauges"
    nwann = size(gauges[1], 2)
    nbvecs = length(overlaps[1])

    wb = kstencil.bweights
    recip_lattice = reciprocal_lattice(kstencil)
    kpoints = kstencil.kpoints

    # Aᵂ can be indexed by Aᵂ[ik][m, n][α], where
    # - ik: kpoint index
    # - m, n: WF indices
    # - α ∈ {1, 2, 3} for x, y, z directions
    Aᵂ = map(1:nkpts) do ik
        Uₖ = gauges[ik]
        Aᵂₖ = zeros(Vec3{eltype(Uₖ)}, nwann, nwann)
        for ib in 1:nbvecs
            ik2 = kstencil.kpb_k[ik][ib]
            Mᴴ = overlaps[ik][ib]
            Uₖ₂ = gauges[ik2]
            Mᵂ = Uₖ' * Mᴴ * Uₖ₂
            G = kstencil.kpb_G[ik][ib]
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
    model::Model, gauges::AbstractVector=model.gauges; kwargs...
)
    return compute_berry_connection_kspace(
        model.kstencil, model.overlaps, gauges; kwargs...
    )
end
