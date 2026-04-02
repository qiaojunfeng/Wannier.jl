
function omega(
    p::CenterSpreadPenalty,
    model::MagModel,
    Uup::AbstractVector{<:AbstractMatrix{T}},
    Udn::AbstractVector{<:AbstractMatrix{T}},
    λs::R,
) where {T<:Complex,R<:Real}
    up = omega(p, model.up, Uup)
    dn = omega(p, model.dn, Udn)
    M = overlap_updn(model, Uup, Udn)
    Ωupdn = omega_updn(M)
    Ωt = up.Ωt + dn.Ωt + λs * Ωupdn
    return SpreadMag(up, dn, Ωupdn, Ωt, M, λs)
end

function omega(p::CenterSpreadPenalty, model::MagModel{T}, λs::T) where {T<:Real}
    return omega(p, model, model.up.gauges, model.dn.gauges, λs)
end

"""
    get_fg!_disentangle(p, model::MagModel, λs)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_disentangle(
    p::AbstractPenalty, model::MagModel{T}, λs::T
) where {T<:Real}
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)
    n_inner = nb * nw + nw^2  # size of XY at each k-point

    function f(XY)
        XY = reshape(XY, (2 * n_inner, nk))  # *2 for spin up and down
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        Ωup = omega(p, model.up.kstencil, model.up.overlaps, Xup, Yup).Ωt
        Ωdn = omega(p, model.dn.kstencil, model.dn.overlaps, Xdn, Ydn).Ωt
        if λs == 0
            Ωupdn = 0
        else
            Ωupdn = omega_updn(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
        end
        return Ωup + Ωdn + λs * Ωupdn
    end

    """size(G) == size(XY)"""
    function g!(G, XY)
        XY = reshape(XY, (2 * n_inner, nk))  # *2 for spin up and down
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        GXup, GYup = omega_grad(
            p, model.up.kstencil, model.up.overlaps, Xup, Yup, model.up.frozen_bands
        )
        GXdn, GYdn = omega_grad(
            p, model.dn.kstencil, model.dn.overlaps, Xdn, Ydn, model.dn.frozen_bands
        )

        # gradient of ↑↓ overlap term
        if λs != 0
            GOXup, GOYup, GOXdn, GOYdn = omega_updn_grad(model, Xup, Yup, Xdn, Ydn)
            GXup += λs * GOXup
            GYup += λs * GOYup
            GXdn += λs * GOXdn
            GYdn += λs * GOYdn
        end

        n = nw^2

        for ik in 1:nk
            G[1:n, ik] = vec(GXup[ik])
            G[(n + 1):n_inner, ik] = vec(GYup[ik])
            G[(n_inner + 1):(n_inner + n), ik] = vec(GXdn[ik])
            G[(n_inner + n + 1):end, ik] = vec(GYdn[ik])
        end

        return nothing
    end

    return f, g!
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
    p::AbstractPenalty,
    model::MagModel{T},
    λs::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
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
    f, g! = get_fg!_disentangle(p, model, λs)

    @info "Initial spread"
    Ω = omega(p, model, λs)
    show(Ω)
    println("\n")

    @info "Initial spread (with states freezed)"
    Ω = omega(p, model, X_Y_to_U(Xup0, Yup0), X_Y_to_U(Xdn0, Ydn0), λs)
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

    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nb, nw)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nb, nw)
    Uupmin = X_Y_to_U(Xupmin, Yupmin)
    Udnmin = X_Y_to_U(Xdnmin, Ydnmin)

    @info "Final spread"
    Ω = omega(p, model, Uupmin, Udnmin, λs)
    show(Ω)

    return Uupmin, Udnmin
end
