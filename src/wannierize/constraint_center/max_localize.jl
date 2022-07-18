using Optim: Optim

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

function pprint(Ω::SpreadCenter)
    println("  WF     center [rx, ry, rz]/Å              spread/Å²  ω  ωc  ωt")

    n_wann = length(Ω.ω)

    for i in 1:n_wann
        @printf(
            "%4d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n",
            i,
            Ω.r[:, i]...,
            Ω.ω[i],
            Ω.ωc[i],
            Ω.ωt[i]
        )
    end

    @printf("Sum spread: Ωt = Ω + Ωc, Ω = ΩI + Ω̃, Ω̃ = ΩOD + ΩD\n")
    @printf("   Ωt  = %11.5f\n", Ω.Ωt)
    @printf("   Ωc  = %11.5f\n", Ω.Ωc)
    @printf("   Ω   = %11.5f\n", Ω.Ω)
    @printf("   ΩI  = %11.5f\n", Ω.ΩI)
    @printf("   ΩOD = %11.5f\n", Ω.ΩOD)
    @printf("   ΩD  = %11.5f\n", Ω.ΩD)
    @printf("   Ω̃   = %11.5f\n", Ω.Ω̃)

    println()

    return nothing
end

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

            MAᵏᵇ .= overlap(M, kpb_k, ik, ikpb) * A[:, :, ikpb]
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
Maximally localize spread functional w.r.t. all kpoints.

On a unitary matrix manifold.
"""
function max_localize_center(
    model::Model{T},
    r₀::Matrix{T},
    λ::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=20,
) where {T<:Real}
    model.n_bands != model.n_wann &&
        error("n_bands != n_wann, run instead disentanglement?")
    size(r₀) != (3, model.n_wann) && error("size(r₀) != (3, n_wann)")

    f, g! = get_fg!_center_maxloc(model, r₀, λ)

    Ωⁱ = omega_center(model.bvectors, model.M, model.A, r₀, λ)
    @info "Initial spread"
    pprint(Ωⁱ)

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
    pprint(Ωᶠ)

    return Amin
end
