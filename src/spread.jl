using LinearAlgebra

export omega, omega_grad, center

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

    # WF center, Cartesian! coordinates, unit Å, 3 * n_wann
    r::Vector{Vec3{T}}

    # frozen_weight::T
    # fix_centers :: Array{Float64,2} #3 x nwannier
end

"""
    omega(bvectors, M, U)

Compute WF spread.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
@views function omega(
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    r = zeros(Vec3{FT}, n_wann)
    r² = zeros(FT, n_wann)

    ΩI::FT = 0.0
    ΩOD::FT = 0.0
    ΩD::FT = 0.0

    b = zeros(FT, 3)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MUᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik in 1:n_kpts
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])

        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]

            MUᵏᵇ .= M[ik][:, :, ib] * U[ikpb]
            # compute-intensive, but should be true
            # ibm = index_bvector(bvectors, ikpb, ik, -kpb_b[:, ib, ik])
            # @assert M[:, :, ib, ik]' ≈ M[:, :, ibm, ikpb]
            Nᵏᵇ .= U[ik]' * MUᵏᵇ
            b .= recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])

            wᵇ = wb[ib]

            ΩI += wᵇ * (n_wann - sum(abs2, Nᵏᵇ))
            ΩOD += wᵇ * sum(abs2, Nᵏᵇ .- diagm(0 => diag(Nᵏᵇ)))

            for n in 1:n_wann
                imlogN = imaglog(Nᵏᵇ[n, n])

                r[n] = r[n] - wᵇ * imlogN * b

                r²[n] += wᵇ * (1 - abs2(Nᵏᵇ[n, n]) + imlogN^2)
                # r²[n] += wᵇ * 2*(1 - real(Nᵏᵇ[n,n]))
            end
        end
    end

    r /= n_kpts
    r² /= n_kpts
    ΩI /= n_kpts
    ΩOD /= n_kpts
    # w_froz /= n_kpts

    # ΩD requires r, so we need different loops
    # However, since ΩD = Ω - ΩI - ΩOD, we can skip these loops
    # for ik in 1:n_kpts
    #     for ib in 1:n_bvecs
    #         ikpb = kpb_k[ib, ik]
    #         Nᵏᵇ .= U[:, :, ik]' * M[:, :, ib, ik] * U[:, :, ikpb]
    #         b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])
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
    ω = r² - map(x -> sum(abs.(x.^2)), r)
    # total Ω
    Ω = sum(ω)
    # Ω += w_froz
    Ω̃ = Ω - ΩI
    ΩD = Ω̃ - ΩOD

    return Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r)
    # return Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r, w_froz)
end

"""
    omega(bvectors, M, U)

Compute WF spread.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
@views function omega(
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Array{Complex{FT}, 3}
) where {FT<:Real}
    n_bands, n_wann, n_kpts = size(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    r = zeros(Vec3{FT}, n_wann)
    r² = zeros(FT, n_wann)

    ΩI::FT = 0.0
    ΩOD::FT = 0.0
    ΩD::FT = 0.0

    b = zeros(FT, 3)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MUᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik in 1:n_kpts
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])

        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]

            MUᵏᵇ .= M[ik][:, :, ib] * U[:, :,ikpb]
            # compute-intensive, but should be true
            # ibm = index_bvector(bvectors, ikpb, ik, -kpb_b[:, ib, ik])
            # @assert M[:, :, ib, ik]' ≈ M[:, :, ibm, ikpb]
            Nᵏᵇ .= U[:, :, ik]' * MUᵏᵇ
            b .= recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])

            wᵇ = wb[ib]

            ΩI += wᵇ * (n_wann - sum(abs2, Nᵏᵇ))
            ΩOD += wᵇ * sum(abs2, Nᵏᵇ .- diagm(0 => diag(Nᵏᵇ)))

            for n in 1:n_wann
                imlogN = imaglog(Nᵏᵇ[n, n])

                r[n] = r[n] - wᵇ * imlogN * b

                r²[n] += wᵇ * (1 - abs2(Nᵏᵇ[n, n]) + imlogN^2)
                # r²[n] += wᵇ * 2*(1 - real(Nᵏᵇ[n,n]))
            end
        end
    end

    r /= n_kpts
    r² /= n_kpts
    ΩI /= n_kpts
    ΩOD /= n_kpts
    # w_froz /= n_kpts

    # ΩD requires r, so we need different loops
    # However, since ΩD = Ω - ΩI - ΩOD, we can skip these loops
    # for ik in 1:n_kpts
    #     for ib in 1:n_bvecs
    #         ikpb = kpb_k[ib, ik]
    #         Nᵏᵇ .= U[:, :, ik]' * M[:, :, ib, ik] * U[:, :, ikpb]
    #         b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])
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
    ω = r² - map(x -> sum(abs.(x.^2)), r)
    # total Ω
    Ω = sum(ω)
    # Ω += w_froz
    Ω̃ = Ω - ΩI
    ΩD = Ω̃ - ΩOD

    return Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r)
    # return Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r, w_froz)
end

"""
    omega(model)

Compute WF spread for `Model`.
"""
omega(model::Model) = omega(model.bvectors, model.M, model.U)

"""
    omega(model, U)

Compute WF spread for `Model` using given `U` gauge.
"""
function omega(model::Model, U::Vector{<:AbstractMatrix})
    return omega(model.bvectors, model.M, U)
end

function Base.show(io::IO, Ω::Spread)
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
function omega_grad(
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Vector{Matrix{Complex{FT}}}, r::Vector{Vec3{FT}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    G = zeros(Complex{FT}, n_bands, n_wann, n_kpts)

    R = zeros(Complex{FT}, n_bands, n_wann)
    T = zeros(Complex{FT}, n_bands, n_wann)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MUᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    q = zeros(FT, n_wann)
    
    @inbounds @views for ik in 1:n_kpts
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])
        # G[1:n_froz, :, ik] = -2 * μ * U[1:n_froz, :, ik]
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]

            mul!(MUᵏᵇ, M[ik][:, :, ib], U[ikpb])
            
            mul!(Nᵏᵇ, U[ik]', MUᵏᵇ)
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])
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

            q .= imaglog.(diag(Nᵏᵇ))
            for ir in 1:n_wann
                q[ir] += r[ir] ⋅ b
            end

            for n in 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                nn = Nᵏᵇ[n, n]
                cnn = conj(nn)
                if abs(nn) < 1e-10
                    display(Nᵏᵇ)
                    println()
                    error("Nᵏᵇ too small! $ik -> $ikpb")
                end

                t = -im * q[n] / nn

                for m in 1:n_bands
                    R[m, n] = -MUᵏᵇ[m, n] * cnn
                    # T[m, n] = -im * MUᵏᵇ[m, n] / (Nᵏᵇ[n, n]) * q[n]
                    T[m, n] = t * MUᵏᵇ[m, n]
                end
            end

            G[:, :, ik] .+= 4 .* wᵇ .* (R .+ T)
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
function omega_grad(
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Array{Complex{FT},3}, r::Vector{Vec3{FT}}
) where {FT<:Real}
    n_bands, n_wann, n_kpts = size(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    G = zeros(Complex{FT}, n_bands, n_wann, n_kpts)

    R = zeros(Complex{FT}, n_bands, n_wann)
    T = zeros(Complex{FT}, n_bands, n_wann)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MUᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    q = zeros(FT, n_wann)
    
    @inbounds @views for ik in 1:n_kpts
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])
        # G[1:n_froz, :, ik] = -2 * μ * U[1:n_froz, :, ik]
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]

            mul!(MUᵏᵇ, M[ik][:, :, ib], U[:, :, ikpb])
            
            mul!(Nᵏᵇ, U[:, :, ik]', MUᵏᵇ)
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])
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

            q .= imaglog.(diag(Nᵏᵇ))
            for ir in 1:n_wann
                q[ir] += r[ir] ⋅ b
            end

            for n in 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                nn = Nᵏᵇ[n, n]
                cnn = conj(nn)
                if abs(nn) < 1e-10
                    display(Nᵏᵇ)
                    println()
                    error("Nᵏᵇ too small! $ik -> $ikpb")
                end

                t = -im * q[n] / nn

                for m in 1:n_bands
                    R[m, n] = -MUᵏᵇ[m, n] * cnn
                    # T[m, n] = -im * MUᵏᵇ[m, n] / (Nᵏᵇ[n, n]) * q[n]
                    T[m, n] = t * MUᵏᵇ[m, n]
                end
            end

            G[:, :, ik] .+= 4 .* wᵇ .* (R .+ T)
        end
    end

    G /= n_kpts

    return G
end

"""
    omega_grad(bvectors, M, U)

Compute gradient of WF spread.

Size of output `dΩ/dU` = `n_bands * n_wann * n_kpts`.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
function omega_grad(
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U
) where {FT<:Real}
    # r = omega(bvectors, M, U).r
    r = center(bvectors, M, U)
    return omega_grad(bvectors, M, U, r)
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
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    wb = bvectors.weights

    loc = zeros(FT, n_kpts)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            Nᵏᵇ .= U[ik]' * M[ik][:, :, ib] * U[ikpb]

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
function center(
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    r = zeros(Vec3{FT}, n_wann)
    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    cache = zeros(Complex{FT}, n_bands, n_wann)
    rt = collect(r')
    # M_ = map(ik -> map(ib -> , 1:n_bvecs), 1:n_kpts)

    @inbounds @views for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            
            mul!(cache, M[ik][:, :, ib], U[ikpb])
            mul!(Nᵏᵇ, U[ik]', cache)
            
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])

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
    center(bvectors, M, U)

Compute WF center in reciprocal space.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
function center(
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Array{Complex{FT},3}
) where {FT<:Real}
    n_bands, n_wann, n_kpts= size(U)
    n_bvecs = size(M[1], 3)
    
    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    r = zeros(Vec3{FT}, n_wann)
    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    cache = zeros(Complex{FT}, n_bands, n_wann)
    rt = collect(r')
    # M_ = map(ik -> map(ib -> , 1:n_bvecs), 1:n_kpts)

    @inbounds @views for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            mul!(cache, M[ik][:, :, ib], U[:,:,ikpb])
            mul!(Nᵏᵇ, U[:,:,ik]', cache)
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])

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
center(model::Model) = center(model.bvectors, model.M, model.U)

"""
    center(model, U)

Compute WF center in reciprocal space for `Model` with given `U` gauge.

# Arguments
- `model`: the `Model`
- `U`: `n_wann * n_wann * n_kpts` array
"""
function center(model::Model, U::Vector{Matrix{T}}) where {T<:Number}
    return center(model.bvectors, model.M, U)
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
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # along x, y, z directions
    R = zeros(Vec3{Complex{FT}}, n_wann, n_wann)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]

            Nᵏᵇ .= U[ik]' * M[ik][:, :, ib] * U[ikpb]
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])

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
position_op(model::Model) = position_op(model.bvectors, model.M, model.U)

"""
    position_op(model, U)

Compute WF postion operator matrix in reciprocal space for `Model` with given `U` gauge.

# Arguments
- `model`: the `Model`
- `U`: `n_wann * n_wann * n_kpts` array
"""
function position_op(model::Model, U::Vector{Matrix{T}}) where {T<:Number}
    return position_op(model.bvectors, model.M, U)
end

"""
    berry_connection(bvectors, M, U)

Compute Berry connection at each kpoint.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
"""
@views function berry_connection(
    bvectors::BVectors{FT}, M::Vector{Matrix{Complex{FT}}}, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # along x, y, z directions
    A = [zeros(Vec3{Complex{FT}}, n_wann, n_wann) for i = 1:n_kpts]
    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]

            Nᵏᵇ .= U[ik]' * M[ik][:, :, ib] * U[ikpb]
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])
            wᵇ = wb[ib]

            for m in 1:n_wann
                for n in 1:n_wann
                    A[ik][m, n] += wᵇ * Nᵏᵇ[m, n] * b

                    if m == n
                        A[ik][m, n] -= wᵇ * b
                    end
                end
            end
        end
    end

    A *= im

    return A
end
