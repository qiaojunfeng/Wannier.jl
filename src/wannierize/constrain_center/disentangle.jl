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
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
    A = X_Y_to_A(X, Y)
    return omega_center(bvectors, M, A, r₀, λ)
end

"""
    omega_center(mode, X, Y, r₀, λ)

Compute WF spread with center penalty, in the `(X, Y)` layout.

# Arguments
- `model`: `Model`
- `X`: `n_wann * n_wann * n_kpts` array
- `Y`: `n_bands * n_wann * n_kpts` array
- `r₀`: `3 * n_wann`, WF centers in cartesian coordinates
- `λ`: penalty strength
"""
function omega_center(
    model::Model{FT},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
    return omega_center(model.bvectors, model.M, X, Y, r₀, λ)
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
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
    frozen::BitMatrix,
    r₀::Matrix{FT},
    λ::FT,
) where {FT<:Real}
    n_kpts = size(Y, 3)

    A = X_Y_to_A(X, Y)
    G = omega_center_grad(bvectors, M, A, r₀, λ)

    GX = zero(X)
    GY = zero(Y)

    for ik in 1:n_kpts
        idx_f = frozen[:, ik]
        n_froz = count(idx_f)

        GX[:, :, ik] = Y[:, :, ik]' * G[:, :, ik]
        GY[:, :, ik] = G[:, :, ik] * X[:, :, ik]'

        GY[idx_f, :, ik] .= 0
        GY[:, 1:n_froz, ik] .= 0
    end

    return GX, GY
end

"""
    get_fg!_center_disentangle(model, r₀, λ=1.0)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_center_disentangle(
    model::Model{T}, r₀::Matrix{T}, λ::T=1.0
) where {T<:Real}
    function f(XY)
        X, Y = XY_to_X_Y(XY, model.n_bands, model.n_wann)
        return omega_center(model.bvectors, model.M, X, Y, r₀, λ).Ωt
    end

    """size(G) == size(XY)"""
    function g!(G, XY)
        X, Y = XY_to_X_Y(XY, model.n_bands, model.n_wann)
        GX, GY = omega_center_grad(model.bvectors, model.M, X, Y, model.frozen_bands, r₀, λ)

        n = model.n_wann^2

        for ik in 1:(model.n_kpts)
            G[1:n, ik] = vec(GX[:, :, ik])
            G[(n + 1):end, ik] = vec(GY[:, :, ik])
        end

        return nothing
    end

    return f, g!
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
    r₀::Matrix{T},
    λ::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
    n_bands = model.n_bands
    n_wann = model.n_wann
    n_kpts = model.n_kpts

    X0, Y0 = A_to_X_Y(model.A, model.frozen_bands)

    # compact storage
    XY0 = X_Y_to_XY(X0, Y0)

    # We have three storage formats:
    # (X, Y): n_wann * n_wann * n_kpts, n_bands * n_wann * n_kpts
    # A: n_bands * n_wann * n_kpts
    # XY: (n_wann * n_wann + n_bands * n_wann) * n_kpts
    f, g! = get_fg!_center_disentangle(model, r₀, λ)

    Ωⁱ = omega_center(model, r₀, λ)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    Ωⁱ = omega_center(model, X0, Y0, r₀, λ)
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

    opt = Optim.optimize(
        f,
        g!,
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
    Amin = X_Y_to_A(Xmin, Ymin)

    Ωᶠ = omega_center(model, Amin, r₀, λ)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Amin
end