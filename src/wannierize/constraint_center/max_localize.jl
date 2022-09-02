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
    r::Matrix{T}

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
            Ω.r[:, i]...,
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
    omega_center(bvectors, M, A, r₀, λ)

Compute WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `A`: `n_wann * n_wann * n_kpts` array
- `r₀`: `3 * n_wann`, WF centers in cartesian coordinates
- `λ`: penalty strength
"""
@views function omega_center(
    bvectors::BVectors{T},
    M::Array{Complex{T},4},
    A::Array{Complex{T},3},
    r₀::Matrix{T},
    λ::T,
) where {T<:Real}
    Ω = omega(bvectors, M, A)
    ωc = λ * dropdims(sum(abs2, Ω.r - r₀; dims=1); dims=1)
    ωt = Ω.ω + ωc
    Ωc = sum(ωc)
    Ωt = Ω.Ω + Ωc
    return SpreadCenter(Ω.Ω, Ω.ΩI, Ω.ΩOD, Ω.ΩD, Ω.Ω̃, Ω.ω, Ω.r, Ωc, Ωt, ωc, ωt)
end

"""
    omega_center_grad(bvectors, M, A, r, r₀, λ)

Compute gradient of WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `A`: `n_wann * n_wann * n_kpts` array
- `r`: `3 * n_wann`, the current WF centers in cartesian coordinates
- `r₀`: `3 * n_wann`, the target WF centers in cartesian coordinates
- `λ`: penalty strength
"""
@views function omega_center_grad(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    r::Matrix{FT},
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
    n_bands, n_wann, n_kpts = size(A)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    G = zeros(Complex{FT}, n_bands, n_wann, n_kpts)

    R = zeros(Complex{FT}, n_bands, n_wann)
    T = zeros(Complex{FT}, n_bands, n_wann)

    b = zeros(FT, 3)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MAᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik in 1:n_kpts
        for ib in 1:n_bvecs
            ikpb = kpb_k[ib, ik]

            MAᵏᵇ .= M[:, :, ib, ik] * A[:, :, ikpb]
            Nᵏᵇ .= A[:, :, ik]' * MAᵏᵇ
            b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])
            wᵇ = wb[ib]

            q = imaglog.(diag(Nᵏᵇ))
            q += r' * b
            # center constraint
            q -= λ * (r - r₀)' * b

            for n in 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                if abs(Nᵏᵇ[n, n]) < 1e-10
                    display(Nᵏᵇ)
                    println()
                    error("Nᵏᵇ too small! $ik -> $ikpb")
                end

                t = -im * q[n] / Nᵏᵇ[n, n]

                for m in 1:n_bands
                    R[m, n] = -MAᵏᵇ[m, n] * conj(Nᵏᵇ[n, n])
                    T[m, n] = t * MAᵏᵇ[m, n]
                end
            end

            G[:, :, ik] .+= 4 * wᵇ .* (R .+ T)
        end
    end

    G /= n_kpts

    return G
end

"""
    omega_center_grad(bvectors, M, A, r₀, λ)

Compute gradient of WF spread with center penalty, for maximal localization.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `A`: `n_wann * n_wann * n_kpts` array
- `r₀`: `3 * n_wann`, the target WF centers in cartesian coordinates
- `λ`: penalty strength
"""
function omega_center_grad(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
    r = center(bvectors, M, A)
    return omega_center_grad(bvectors, M, A, r, r₀, λ)
end

"""
    get_fg!_center_maxloc(model, r₀, λ=1.0)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_center_maxloc(model::Model{T}, r₀::Matrix{T}, λ::T=1.0) where {T<:Real}
    f(A) = omega_center(model.bvectors, model.M, A, r₀, λ).Ωt

    function g!(G, A)
        r = center(model.bvectors, model.M, A)
        G .= omega_center_grad(model.bvectors, model.M, A, r, r₀, λ)
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
    r₀::Matrix{T},
    λ::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
    model.n_bands != model.n_wann &&
        error("n_bands != n_wann, run instead disentanglement?")
    size(r₀) != (3, model.n_wann) && error("size(r₀) != (3, n_wann)")

    f, g! = get_fg!_center_maxloc(model, r₀, λ)

    Ωⁱ = omega_center(model.bvectors, model.M, model.A, r₀, λ)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    kManif = Optim.Stiefel_SVD()
    Manif = Optim.PowerManifold(kManif, (model.n_wann, model.n_wann), (model.n_kpts,))

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    Ainit = deepcopy(model.A)

    opt = Optim.optimize(
        f,
        g!,
        Ainit,
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

    Amin = Optim.minimizer(opt)

    Ωᶠ = omega_center(model.bvectors, model.M, Amin, r₀, λ)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Amin
end
