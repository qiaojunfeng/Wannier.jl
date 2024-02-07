
function omega(
    model::MagModel,
    Uup::AbstractVector,
    Udn::AbstractVector,
    r₀::AbstractVector,
    λc::Real,
    λs::Real,
)
    up = omega_center(model.up, Uup; r₀, λc)
    dn = omega_center(model.dn, Udn; r₀, λc)
    M = overlap_updn(model, Uup, Udn)
    Ωupdn = omega_updn(M)
    Ωt = up.Ω + dn.Ω + λs * Ωupdn
    return SpreadMag(up, dn, Ωupdn, Ωt, M, λs)
end

function omega(model::MagModel, r₀::AbstractVector, λc::Real, λs::Real)
    return omega(model, model.up.gauges, model.dn.gauges, r₀, λc, λs)
end

"""
    get_fg!_center_disentangle(model::MagModel, r₀, λc, λs)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_center_disentangle(
    model::MagModel{T}, r₀::Vector{Vec3{T}}, λc::T, λs::T
) where {T<:Real}
    nbands = n_bands(model.up)
    nwann = n_wannier(model.up)
    nkpts = n_kpoints(model.up)
    ninner = nbands * nwann + nwann^2  # size of XY at each k-point

    function f(XY)
        XY = reshape(XY, (2 * ninner, nkpts))  # *2 for spin up and down
        XYup = @view XY[1:ninner, :]
        XYdn = @view XY[(ninner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nbands, nwann)
        Xdn, Ydn = XY_to_X_Y(XYdn, nbands, nwann)
        Ωup = omega_center(model.up.kstencil, model.up.overlaps, Xup, Yup, r₀, λc).Ω
        Ωdn = omega_center(model.dn.kstencil, model.dn.overlaps, Xdn, Ydn, r₀, λc).Ω
        if λs == 0
            Ωupdn = 0
        else
            Ωupdn = omega_updn(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
        end
        return Ωup + Ωdn + λs * Ωupdn
    end

    """size(G) == size(XY)"""
    function g!(G, XY)
        XY = reshape(XY, (2 * ninner, nkpts))  # *2 for spin up and down
        XYup = @view XY[1:ninner, :]
        XYdn = @view XY[(ninner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nbands, nwann)
        Xdn, Ydn = XY_to_X_Y(XYdn, nbands, nwann)
        GXup, GYup = omega_center_grad(
            model.up.kstencil, model.up.overlaps, Xup, Yup, model.up.frozen_bands, r₀, λc
        )
        GXdn, GYdn = omega_center_grad(
            model.dn.kstencil, model.dn.overlaps, Xdn, Ydn, model.dn.frozen_bands, r₀, λc
        )

        # gradient of ↑↓ overlap term
        if λs != 0
            GOXup, GOYup, GOXdn, GOYdn = omega_updn_grad(model, Xup, Yup, Xdn, Ydn)
            GXup += λs * GOXup
            GYup += λs * GOYup
            GXdn += λs * GOXdn
            GYdn += λs * GOYdn
        end

        n = nwann^2

        for ik in 1:nkpts
            G[1:n, ik] = vec(GXup[:, :, ik])
            G[(n + 1):ninner, ik] = vec(GYup[:, :, ik])
            G[(ninner + 1):(ninner + n), ik] = vec(GXdn[:, :, ik])
            G[(ninner + n + 1):end, ik] = vec(GYdn[:, :, ik])
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
- `r₀`: `3 * n_wann`, WF centers in cartesian coordinates
- `λc`: Lagrange multiplier of the center term
- `λs`: Lagrange multiplier of the spin-up and spin-down overlap term

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function disentangle_center(
    model::MagModel{T},
    r₀::Vector{Vec3{T}},
    λc::T=1.0,
    λs::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
    nbands = n_bands(model.up)
    nwann = n_wannier(model.up)
    nkpts = n_kpoints(model.up)

    @assert n_bands(model.dn) == nbands
    @assert n_wannier(model.dn) == nwann
    @assert n_kpoints(model.dn) == nkpts

    XYk_up_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nwann, nwann), (nbands, nwann)
    )
    XYk_dn_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nwann, nwann), (nbands, nwann)
    )
    n_inner = nwann^2 + nbands * nwann
    XYkManif = Optim.ProductManifold(XYk_up_Manif, XYk_dn_Manif, (n_inner,), (n_inner,))
    XYManif = Optim.PowerManifold(XYkManif, (2 * n_inner,), (nkpts,))

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
    f, g! = get_fg!_center_disentangle(model, r₀, λc, λs)

    Ω = omega(model, r₀, λc, λs)
    @info "Initial spread" Ω
    println("\n")

    Ω = omega(model, X_Y_to_U(Xup0, Yup0), X_Y_to_U(Xdn0, Ydn0), r₀, λc, λs)
    @info "Initial spread (with states freezed)" Ω
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
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nbands, nwann)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nbands, nwann)
    Uupmin = X_Y_to_U(Xupmin, Yupmin)
    Udnmin = X_Y_to_U(Xdnmin, Ydnmin)

    @info "Final spread"
    Ω = omega(model, Uupmin, Udnmin, r₀, λc, λs)
    show(Ω)

    return Uupmin, Udnmin
end
