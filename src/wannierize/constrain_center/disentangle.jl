using LinearAlgebra
using Optim: Optim

export disentangle_center

"""
    omega_center(bvectors, M, X, Y, r₀, λ)

Compute WF spread with center penalty, in the `(X, Y)` layout.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `X`: `n_wann * n_wann * n_kpts` array
- `Y`: `n_bands * n_wann * n_kpts` array
- `r₀`: `3 * n_wann`, WF centers in cartesian coordinates
- `λ`: penalty strength
"""
function omega_center(
    bvectors::BVectors{FT},
    M::Vector{Array{Complex{FT},3}},
    X::Vector{Matrix{Complex{FT}}},
    Y::Vector{Matrix{Complex{FT}}},
    r₀::Vector{Vec3{FT}},
    λ::FT,
) where {FT<:Real}
    U = X_Y_to_U(X, Y)
    return omega_center(bvectors, M, U; r₀, λ)
end

"""
    omega_center_grad(bvectors, M, X, Y, frozen, r₀, λ)

Compute gradient of WF spread with center penalty, in the `(X, Y)` layout.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `X`: `n_wann * n_wann * n_kpts` array
- `Y`: `n_bands * n_wann * n_kpts` array
- `frozen`: `n_bands * n_kpts` array for frozen bands
- `r₀`: `3 * n_wann`, WF centers in cartesian coordinates
- `λ`: penalty strength
"""
function omega_center_grad(
    bvectors::BVectors{FT},
    M::Vector{Array{Complex{FT},3}},
    X::Vector{Matrix{Complex{FT}}},
    Y::Vector{Matrix{Complex{FT}}},
    frozen::Vector{BitVector},
    r₀::Vector{Vec3{FT}},
    λ::FT,
) where {FT<:Real}
    return omega_grad(center_penalty(r₀, λ), bvectors, M, X, Y, frozen)
end

"""
    get_fg!_center_disentangle(model, r₀, λ=1.0)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_center_disentangle(
    model::Model{T}, r₀::Vector{Vec3{T}}, λ::T=1.0
) where {T<:Real}
    cache = Cache(model)
    
    function fg!(Ω, G, XY)
        
        X, Y = XY_to_X_Y!(cache.X, cache.Y, XY)
        U = X_Y_to_U!(cache.U, X, Y)
        compute_MUᵏᵇ_Nᵏᵇ!(cache, model.bvectors, model.M, U)

        if G !== nothing
            # TODO Optimize this!
            G_ = omega_center_grad!(cache, model.bvectors, model.M; r₀, λ)
            GX, GY = GU_to_GX_GY(G_, X, Y, model.frozen_bands)

            n = model.n_wann^2

            @inbounds for ik in 1:(model.n_kpts)
                for i in eachindex(GX[ik])
                    G[i, ik] = GX[ik][i]
                end
                for i in eachindex(GY[ik])
                    G[n + i, ik] = GY[ik][i]
                end
            end
        end
        if Ω !== nothing
            return omega_center(omega!(cache, model.bvectors, model.M); r₀, λ).Ωt
        end
    end
    return fg!
end

"""
    disentangle(model, r₀, λ=1.0; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Run disentangle on the `Model` with center penalty.

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
function disentangle_center(
    model::Model{T},
    r₀::Vector{Vec3{T}},
    λ::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
    n_bands = model.n_bands
    n_wann = model.n_wann
    n_kpts = model.n_kpts

    X0, Y0 = U_to_X_Y(model.U, model.frozen_bands)

    # compact storage
    XY0 = X_Y_to_XY(X0, Y0)
    # We have three storage formats:
    # (X, Y): n_wann * n_wann * n_kpts, n_bands * n_wann * n_kpts
    # U: n_bands * n_wann * n_kpts
    # XY: (n_wann * n_wann + n_bands * n_wann) * n_kpts
    fg! = get_fg!_center_disentangle(model, r₀, λ)

    Ωⁱ = omega_center(model; r₀, λ)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    Ωⁱ = omega_center(model, X_Y_to_U(X0, Y0); r₀, λ)
    @info "Initial spread (with states freezed)"
    show(Ωⁱ)
    println("\n")

    # need QR orthogonalization rather than SVD to preserve the sparsity structure of Y
    XYkManif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (n_wann, n_wann), (n_bands, n_wann)
    )
    XYManif = Optim.PowerManifold(XYkManif, (n_wann^2 + n_bands * n_wann,), (n_kpts,))

    # stepsize_mult = 1
    # step = 0.5/(4*8*p.wb)*(p.N1*p.N2*p.N3)*stepsize_mult
    # ls = LineSearches.Static(step)
    ls = Optim.HagerZhang()
    # ls = LineSearches.BackTracking()

    # meth = Optim.GradientDescent
    # meth = Optim.ConjugateGradient
    meth = Optim.LBFGS

    opt = Optim.optimize(Optim.only_fg!(fg!),
        XY0,
        meth(; manifold=XYManif, linesearch=ls, m=history_size),
        Optim.Options(;
            show_trace=true,
            iterations=max_iter,
            f_tol=f_tol,
            g_tol=g_tol,
            allow_f_increases=true,
        ),
    )
    display(opt)

    XYmin = Optim.minimizer(opt)

    Xmin, Ymin = XY_to_X_Y(XYmin, n_bands, n_wann)
    Umin = X_Y_to_U(Xmin, Ymin)

    Ωᶠ = omega_center(model, Umin; r₀, λ)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Umin
end
