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

struct DisentangleCache{T}
    X::Vector{Matrix{Complex{T}}}
    Y::Vector{Matrix{Complex{T}}}
    U::Vector{Matrix{Complex{T}}}
    G::Array{Complex{T}, 3}
    r::Vector{Vec3{T}}
    Nᵏᵇ::Vector{Vector{Matrix{Complex{T}}}}
    MUᵏᵇ::Vector{Vector{Matrix{Complex{T}}}}
end

function DisentangleCache(bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Vector{Matrix{Complex{FT}}}) where {FT}
    nk = length(U)
    nb, nw = size(U[1])
    nbvecs = size(M[1], 3)
    
    X = [zeros(Complex{FT}, nw, nw) for i = 1:nk]
    Y = [zeros(Complex{FT}, nb, nw) for i = 1:nk]
    U = [zeros(Complex{FT}, nb, nw) for i = 1:nk]
    G = zeros(Complex{FT}, nb, nw, nk)
    r = zeros(Vec3{FT}, nw)
    
    Nᵏᵇ  = [[zeros(Complex{FT}, nw, nw) for ib=1:nbvecs] for i = 1:nk]
    MUᵏᵇ = [[zeros(Complex{FT}, nb, nw) for ib=1:nbvecs] for i = 1:nk]
    
    return DisentangleCache(X, Y, U, G, r, Nᵏᵇ, MUᵏᵇ)
end

DisentangleCache(model::Model) = DisentangleCache(model.bvectors, model.M, model.U) 

function compute_MUᵏᵇ_Nᵏᵇ!(MUᵏᵇ, Nᵏᵇ, bvectors, M, U)
    kpb_k = bvectors.kpb_k
    n_bvecs = length(kpb_k[1])

    @inbounds for ik in 1:length(U)
        Ut = U[ik]'
        Mk = M[ik]
        MUk = MUᵏᵇ[ik]
        Nk = Nᵏᵇ[ik]
        
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            MUkb = MUk[ib]
            Ukpb = U[ikpb]
            Nkb = Nk[ib]

            mul!(MUkb, view(Mk,:, :, ib), Ukpb)
            
            mul!(Nkb, Ut, MUkb)
        end
    end
    return MUᵏᵇ, Nᵏᵇ
end
compute_MUᵏᵇ_Nᵏᵇ!(cache::DisentangleCache, bvectors, M, U) = compute_MUᵏᵇ_Nᵏᵇ!(cache.MUᵏᵇ, cache.Nᵏᵇ, bvectors, M, U)

function omega!(cache::DisentangleCache,
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Vector{Matrix{Complex{FT}}}
) where {FT<:Real}

    r = cache.r
    fill!(r, zero(eltype(r)))
    
    Nᵏᵇ  = cache.Nᵏᵇ
    MUᵏᵇ = cache.MUᵏᵇ

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

    r² = zeros(FT, n_wann)

    ΩI::FT = 0.0
    ΩOD::FT = 0.0
    ΩD::FT = 0.0

    for ik in 1:n_kpts
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])
        MUk = MUᵏᵇ[ik]
        Nk = Nᵏᵇ[ik]

        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]

            ikpb = kpb_k[ik][ib]
            MUkb = MUk[ib]
            Nkb  = Nk[ib] 
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])

            wᵇ = wb[ib]

            wb_b = wᵇ * b
            
            ts = zero(FT)
            ts2 = zero(FT)
            for i in axes(Nkb, 2)
                for j in axes(Nkb, 1)
                    nt =  Nkb[j, i]
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
            
            ΩI += wᵇ * (n_wann - ts)
            
            ΩOD += wᵇ * ts2
            
        end
    end

    r = map(x -> x./n_kpts, r)
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
    omega(model, [U])
    omega(bvectors, M, U)

Compute WF spread for a [`Model`](@ref), potentially for a given gauge `U`, or by explicitely giving
`bvectors` and `M`.
In case of the first `bvectors = model.bvectors` and `M = model.M`.
"""
omega(model::Model)                              = omega(model, model.U)
omega(model::Model, U::Vector{<:AbstractMatrix}) = omega(model.bvectors, model.M, U)
function omega(bvectors::BVectors, M, U)
    cache = DisentangleCache(bvectors, M, U)
    compute_MUᵏᵇ_Nᵏᵇ!(cache, bvectors, M, U)
    return omega!(cache, bvectors, M, U)
end
function omega(bvectors::BVectors, M, X, Y)
    U = X_Y_to_U(X, Y)
    return omega(bvectors, M, U)
end

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

    r = map(x -> x./n_kpts, r)
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

function omega_grad!(cache::DisentangleCache{T}, bvectors, M, U) where {T}
    # This mutates cache.G and cache.Mkb
    G = cache.G
    fill!(G, 0)
    r = cache.r
    Nᵏᵇ  = cache.Nᵏᵇ
    MUᵏᵇ = cache.MUᵏᵇ
    
    n_bands, n_wann, n_kpts = size(G)

    n_bvecs = size(M[1], 3)

    center!(r, Nᵏᵇ, bvectors)
    
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

    @inbounds for ik in 1:n_kpts
        # w_froz -= μ * sum(abs2, U[1:n_froz, :, ik])
        # G[1:n_froz, :, ik] = -2 * μ * U[1:n_froz, :, ik]
        MUk = MUᵏᵇ[ik]
        Nk = Nᵏᵇ[ik]
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            MUkb = MUk[ib]
            Nkb  = Nk[ib] 
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

            for n in 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                nn = Nkb[n, n]
                
                # TODO: This check can be done somewherhe else adds 12% of time
                # if abs(nn) < 1e-10
                #     display(Nk[ib])
                #     println()
                #     error("Nᵏᵇ too small! $ik -> $ikpb")
                # end

                q = imaglog(nn) + r[n] ⋅ b

                t = -im * q / nn

                cnn = conj(nn)
                for m in 1:n_bands
                    # T[m, n] = -im * MUᵏᵇ[m, n] / (Nᵏᵇ[n, n]) * q[n]
                    MUkb[m, n] *= (t - cnn)
                end
            end

            view(G,:, :, ik) .+= 4 .* wᵇ .* MUk[ib]
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
function omega_grad(bvectors::BVectors, M, U)
    cache = DisentangleCache(bvectors, M, U)
    compute_MUᵏᵇ_Nᵏᵇ!(cache, bvectors, M, U)
    omega_grad!(cache, bvectors, M, U)
end

function omega_grad(bvectors::BVectors, M, X, Y,frozen)
    U = X_Y_to_U(X, Y)
    G = omega_grad(bvectors, M, U)
    return GU_to_GX_GY(G, X, Y, frozen)
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
    bvectors::BVectors{FT}, M::Vector{Array{Complex{FT},3}}, U::Array{Complex{FT},3}
) where {FT<:Real}

    r = center(bvectors, M, U)
    
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
function center(bvectors::BVectors, M, U)
    cache = DisentangleCache(bvectors, M, U)
    compute_MUᵏᵇ_Nᵏᵇ!(cache, bvectors, M, U)
    return center!(cache.r, cache.Nᵏᵇ, bvectors)
end
function center!(r::Vector{<:Vec3}, Nᵏᵇ, bvectors)
    fill!(r, zero(eltype(r)))
    n_wann = length(r)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    @inbounds for (ik, Nk) in enumerate(Nᵏᵇ)
        k = kpoints[ik] 
        for (ib, Nb) in enumerate(Nk)
            ikpb = kpb_k[ik][ib]
            
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - k)

            w = wb[ib]
            
            for n in 1:n_wann
                fac = w * imaglog(Nb[n, n])
                r[n] -= b * fac
            end
        end
    end

    r ./= length(Nᵏᵇ)

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
