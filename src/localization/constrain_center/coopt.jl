function omega(
        terms::Tuple,
        model::MagModel,
        Uup::AbstractVector{<:AbstractMatrix{T}},
        Udn::AbstractVector{<:AbstractMatrix{T}},
        λs::R,
    ) where {T <: Complex, R <: Real}
    center_term = _find_center_term(terms)
    isnothing(center_term) && error("CenterConstraintTerm is required for constrained-center magnetic omega")
    up = omega_center(omega(model.up, Uup); r₀ = center_term.r0, λ = center_term.λ)
    dn = omega_center(omega(model.dn, Udn); r₀ = center_term.r0, λ = center_term.λ)
    M = overlap_updn(model, Uup, Udn)
    Ωupdn = omega_updn(M)
    Ωt = up.Ωt + dn.Ωt + λs * Ωupdn
    return SpreadMag(up, dn, Ωupdn, Ωt, M, λs)
end

function omega(terms::Tuple, model::MagModel{T}, λs::T) where {T <: Real}
    return omega(terms, model, model.up.gauges, model.dn.gauges, λs)
end

"""
    get_fg!_disentangle(p, model::MagModel, λs)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_disentangle(
        terms::Tuple, model::MagModel{T}, λs::T
    ) where {T <: Real}
    problem = LocalizationProblem(terms, model, :mag_disentangle_center; lambda = λs)
    return build_fg!(problem)
end

"""
    disentangle(model; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Run disentangle on a `MagModel`, with center constraints.

# Arguments
- `model`: MagModel
- `λs`: Lagrange multiplier of the spin-up and spin-down overlap term

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function disentangle(
    terms::Tuple,
        model::MagModel{T},
        λs::T = 1.0;
        f_tol::T = 1.0e-7,
        g_tol::T = 1.0e-5,
        max_iter::Int = 200,
        history_size::Int = 3,
    ) where {T <: Real}
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)
    @assert n_bands(model.dn) == nb
    @assert n_wannier(model.dn) == nw
    @assert n_kpoints(model.dn) == nk

    XYk_up_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw)
    )
    XYk_dn_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw)
    )
    n_inner = nw^2 + nb * nw
    XYkManif = Optim.ProductManifold(XYk_up_Manif, XYk_dn_Manif, (n_inner,), (n_inner,))
    XYManif = Optim.PowerManifold(XYkManif, (2 * n_inner,), (nk,))

    Xup0, Yup0 = U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
    # compact storage
    XYup0 = X_Y_to_XY(Xup0, Yup0)
    XYdn0 = X_Y_to_XY(Xdn0, Ydn0)
    XY0 = vcat(XYup0, XYdn0)

    # We have three storage formats:
    # (X, Y): n_wann * n_wann * n_kpts, n_bands * n_wann * n_kpts
    # U: n_bands * n_wann * n_kpts
    # XY: (n_wann * n_wann + n_bands * n_wann) * n_kpts
    f, g! = get_fg!_disentangle(terms, model, λs)

    @info "Initial spread"
    Ω = omega(terms, model, λs)
    show(Ω)
    println("\n")

    @info "Initial spread (with states freezed)"
    Ω = omega(terms, model, X_Y_to_U(Xup0, Yup0), X_Y_to_U(Xdn0, Ydn0), λs)
    show(Ω)
    println("\n")

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
        meth(; manifold = XYManif, linesearch = ls, m = history_size),
        Optim.Options(;
            show_trace = true,
            iterations = max_iter,
            f_tol = f_tol,
            g_tol = g_tol,
            allow_f_increases = true,
        ),
    )
    display(opt)

    XYmin = Optim.minimizer(opt)

    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nb, nw)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nb, nw)
    Uupmin = X_Y_to_U(Xupmin, Yupmin)
    Udnmin = X_Y_to_U(Xdnmin, Ydnmin)

    @info "Final spread"
    Ω = omega(terms, model, Uupmin, Udnmin, λs)
    show(Ω)

    return Uupmin, Udnmin
end

