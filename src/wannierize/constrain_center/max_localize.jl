using Optim: Optim

export max_localize_center

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

function Base.show(io::IO, Ω::SpreadCenter)
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
function omega_center(args...; r₀::Vector{Vec3{T}}, λ::T) where {T<:Real}
    Ω = omega(args...)
    ωc = λ .* map(i -> (t = Ω.r[i] - r₀[i]; sum(t.^2)), 1:length(r₀))
    ωt = Ω.ω + ωc
    Ωc = sum(ωc)
    Ωt = Ω.Ω + Ωc
    return SpreadCenter(Ω.Ω, Ω.ΩI, Ω.ΩOD, Ω.ΩD, Ω.Ω̃, Ω.ω, Ω.r, Ωc, Ωt, ωc, ωt)
end

"""
    omega_center_grad(bvectors, M, U, r, r₀, λ)

Compute gradient of WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
- `r`: `3 * n_wann`, the current WF centers in cartesian coordinates
- `r₀`: `3 * n_wann`, the target WF centers in cartesian coordinates
- `λ`: penalty strength
"""
@views function omega_center_grad(
    bvectors::BVectors{FT},
    M::Vector{Array{Complex{FT},3}},
    U::Array{Complex{FT}, 3},
    r::Vector{Vec3{FT}},
    r₀::Vector{Vec3{FT}},
    λ::FT,
) where {FT<:Real}
    n_bands, n_wann, n_kpts = size(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    G = zeros(Complex{FT}, n_bands, n_wann, n_kpts)

    R = zeros(Complex{FT}, n_bands, n_wann)
    T = zeros(Complex{FT}, n_bands, n_wann)


    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MUᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            #TODO optimize
            MUᵏᵇ .= M[ik][:, :, ib] * U[:,:,ikpb]
            Nᵏᵇ .= U[:,:,ik]' * MUᵏᵇ
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])
            wᵇ = wb[ib]

            q = imaglog.(diag(Nᵏᵇ))
            for ir = 1:n_wann
                q[ir] += r[ir] ⋅ b
            # center constraint
                q[ir] -= λ * (r[ir] - r₀[ir]) ⋅ b
            end

            for n in 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                if abs(Nᵏᵇ[n, n]) < 1e-10
                    display(Nᵏᵇ)
                    println()
                    error("Nᵏᵇ too small! $ik -> $ikpb")
                end

                t = -im * q[n] / Nᵏᵇ[n, n]

                for m in 1:n_bands
                    R[m, n] = -MUᵏᵇ[m, n] * conj(Nᵏᵇ[n, n])
                    T[m, n] = t * MUᵏᵇ[m, n]
                end
            end

            G[:, :, ik] .+= 4 * wᵇ .* (R .+ T)
        end
    end

    G /= n_kpts
    return G
end

"""
    omega_center_grad(bvectors, M, U, r, r₀, λ)

Compute gradient of WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
- `r`: `3 * n_wann`, the current WF centers in cartesian coordinates
- `r₀`: `3 * n_wann`, the target WF centers in cartesian coordinates
- `λ`: penalty strength
"""
@views function omega_center_grad(
    bvectors::BVectors{FT},
    M::Vector{Array{Complex{FT},3}},
    U::Vector{Matrix{Complex{FT}}},
    r::Vector{Vec3{FT}},
    r₀::Vector{Vec3{FT}},
    λ::FT,
) where {FT<:Real}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    n_bvecs = size(M[1], 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    G = [zeros(Complex{FT}, n_bands, n_wann) for i = 1:n_kpts]

    R = zeros(Complex{FT}, n_bands, n_wann)
    T = zeros(Complex{FT}, n_bands, n_wann)


    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MUᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            #TODO optimize
            MUᵏᵇ .= M[ik][:, :, ib] * U[ikpb]
            Nᵏᵇ .= U[ik]' * MUᵏᵇ
            b = recip_lattice * (kpoints[ikpb] + kpb_b[ik][ib] - kpoints[ik])
            wᵇ = wb[ib]

            q = imaglog.(diag(Nᵏᵇ))
            for ir = 1:n_wann
                q[ir] += r[ir] ⋅ b
            # center constraint
                q[ir] -= λ * (r[ir] - r₀[ir]) ⋅ b
            end

            for n in 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                if abs(Nᵏᵇ[n, n]) < 1e-10
                    display(Nᵏᵇ)
                    println()
                    error("Nᵏᵇ too small! $ik -> $ikpb")
                end

                t = -im * q[n] / Nᵏᵇ[n, n]

                for m in 1:n_bands
                    R[m, n] = -MUᵏᵇ[m, n] * conj(Nᵏᵇ[n, n])
                    T[m, n] = t * MUᵏᵇ[m, n]
                end
            end

            G[ik] .+= 4 * wᵇ .* (R .+ T)
        end
    end

    G ./= n_kpts
    return G
end

"""
    omega_center_grad(bvectors, M, U, r₀, λ)

Compute gradient of WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `U`: `n_wann * n_wann * n_kpts` array
- `r₀`: `3 * n_wann`, the target WF centers in cartesian coordinates
- `λ`: penalty strength
"""
function omega_center_grad(
    bvectors::BVectors{FT},
    M::Vector,
    U,
    r₀::Vector,
    λ::FT,
) where {FT<:Real}
    r = center(bvectors, M, U)
    return omega_center_grad(bvectors, M, U, r, r₀, λ)
end

"""
    get_fg!_center_maxloc(model, r₀, λ=1.0)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_center_maxloc(model::Model{T}, r₀::Vector{Vec3{T}}, λ::T=1.0) where {T<:Real}
    f(U) = omega_center(model.bvectors, model.M, U; r₀, λ).Ωt

    function g!(G, U)
        r = center(model.bvectors, model.M, U)
        G .= omega_center_grad(model.bvectors, model.M, U, r, r₀, λ)
        return nothing
    end

    return f, g!
end

"""
    max_localize_center(model, r₀, λ=1.0; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Maximally localize spread functional with center constraint on a unitary matrix manifold.

# Arguments
- `model`: model
- `r₀`: `3 * n_wann`, WF centers in cartesian coordinates
- `λ`: penalty strength

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function max_localize_center(
    model::Model{T},
    r₀::Vector,
    λ::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
    model.n_bands != model.n_wann &&
        error("n_bands != n_wann, run instead disentanglement?")
    length(r₀) !=  model.n_wann && error("length(r₀) !=  n_wann")

    f, g! = get_fg!_center_maxloc(model, r₀, λ)

    Ωⁱ = omega_center(model.bvectors, model.M, model.U; r₀, λ)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    kManif = Optim.Stiefel_SVD()
    Manif = Optim.PowerManifold(kManif, (model.n_wann, model.n_wann), (model.n_kpts,))

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    Uinit = [model.U[ik][ib, ic] for ib=1:size(model.U[1],1), ic = 1:size(model.U[1],2), ik = 1:length(model.U)]

    opt = Optim.optimize(
        f,
        g!,
        Uinit,
        meth(; manifold=Manif, linesearch=ls, m=history_size),
        Optim.Options(;
            show_trace=true,
            iterations=max_iter,
            f_tol=f_tol,
            g_tol=g_tol,
            allow_f_increases=true,
        ),
    )
    display(opt)

    Umin = Optim.minimizer(opt)

    Ωᶠ = omega_center(model.bvectors, model.M, Umin; r₀, λ)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Umin
end
