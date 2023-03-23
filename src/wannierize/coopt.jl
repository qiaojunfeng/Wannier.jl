struct MagModel{T<:Real}
    # submodel for spin up
    up::Model{T}

    # submodel for spin down
    dn::Model{T}

    # < u_nk^↑ | u_mk^↓ >, size: (n_bands, n_bands, n_kpts)
    O::Array{Complex{T},3}
end

# function MagModel(up::Model{T}, dn::Model{T}) where {T<:Real}
#     @assert up.n_bands == dn.n_bands
#     @assert up.n_wann == dn.n_wann
#     @assert up.n_kpts == dn.n_kpts
#     @assert up.n_bvecs == dn.n_bvecs
#     @assert up.kgrid == dn.kgrid
#     @assert up.kpoints == dn.kpoints
#     @assert up.bvectors == dn.bvectors
#     # @assert up.frozen_bands == dn.frozen_bands

#     MagModel(up, dn)
# end

function updn_overlap(
    model::MagModel{T}, Uup::Array{Complex{T},3}, Udn::Array{Complex{T},3}
) where {T<:Real}
    n_bands, n_wann, n_kpts = size(Uup)
    O = zeros(real(eltype(model.O)), n_wann, n_wann)

    for ik in 1:n_kpts
        Uupk = @view Uup[:, :, ik]
        Udnk = @view Udn[:, :, ik]
        Ok = @view model.O[:, :, ik]
        upOkdn = Uupk' * Ok * Udnk
        O += abs2.(upOkdn)
    end
    return O / n_kpts
end

function updn_overlap(model::MagModel)
    return updn_overlap(model, model.up.U, model.dn.U)
end

function updn_overlap_grad(model::MagModel, Xup, Yup, Xdn, Ydn)
    n_bands, n_wann, n_kpts = size(Yup)

    T = eltype(Xup)
    S = zeros(T, n_wann, n_wann)
    GXup = zeros(T, size(Xup))
    GXdn = zeros(T, size(Xdn))
    GYup = zeros(T, size(Yup))
    GYdn = zeros(T, size(Ydn))

    for ik in 1:n_kpts
        Uupk = Yup[:, :, ik] * Xup[:, :, ik]
        Udnk = Ydn[:, :, ik] * Xdn[:, :, ik]
        Ok = @view model.O[:, :, ik]
        upOkdn = Uupk' * Ok * Udnk
        S += upOkdn

        GXup[:, :, ik] = Yup[:, :, ik]' * Ok * Udnk
        GYup[:, :, ik] = Ok * Udnk * Xup[:, :, ik]'
        GXdn[:, :, ik] = Ydn[:, :, ik]' * Ok' * Uupk
        GYdn[:, :, ik] = Ok' * Uupk * Xdn[:, :, ik]'
    end
    for ik in 1:n_kpts
        GXup[:, :, ik] .= GXup[:, :, ik] * S'
        GYup[:, :, ik] .= GYup[:, :, ik] * S'
        GXdn[:, :, ik] .= GXdn[:, :, ik] * S
        GYdn[:, :, ik] .= GYdn[:, :, ik] * S
    end
    return GXup / n_kpts^2, GYup / n_kpts^2, GXdn / n_kpts^2, GYdn / n_kpts^2
end

"""
    get_fg!_disentangle(model::MagModel)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_disentangle(model::MagModel, n_bands, n_wann, n_kpts, λ::Real=1.0)
    function f(XY)
        n_inner = n_bands * n_wann + n_wann^2  # size of XY at each k-point
        XY = reshape(XY, (2 * n_inner, n_kpts))  # *2 for spin up and down
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, n_bands, n_wann)
        Xdn, Ydn = XY_to_X_Y(XYdn, n_bands, n_wann)
        Ωup = omega(model.up.bvectors, model.up.M, Xup, Yup).Ω
        Ωdn = omega(model.dn.bvectors, model.dn.M, Xdn, Ydn).Ω
        Oupdn = updn_overlap(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
        # I am using minus sign here because the optimizer is minimizing total,
        # thus maximizing the ↑↓ overlap.
        return Ωup + Ωdn - λ * sum(diag(Oupdn))
    end

    """size(G) == size(XY)"""
    function g!(G, XY)
        n_inner = n_bands * n_wann + n_wann^2  # size of XY at each k-point
        XY = reshape(XY, (2 * n_inner, n_kpts))  # *2 for spin up and down
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, n_bands, n_wann)
        Xdn, Ydn = XY_to_X_Y(XYdn, n_bands, n_wann)
        GXup, GYup = omega_grad(
            model.up.bvectors, model.up.M, Xup, Yup, model.up.frozen_bands
        )
        GXdn, GYdn = omega_grad(
            model.dn.bvectors, model.dn.M, Xdn, Ydn, model.dn.frozen_bands
        )

        # gradient of ↑↓ overlap term
        GOXup, GOYup, GOXdn, GOYdn = updn_overlap_grad(model, Xup, Yup, Xdn, Ydn)
        GXup -= λ * GOXup
        GYup -= λ * GOYup
        GXdn -= λ * GOXdn
        GYdn -= λ * GOYdn

        n = n_wann^2

        for ik in 1:(model.up.n_kpts)
            G[1:n, ik] = vec(GXup[:, :, ik])
            G[(n + 1):n_inner, ik] = vec(GYup[:, :, ik])
            G[(n_inner + 1):(n_inner + n), ik] = vec(GXdn[:, :, ik])
            G[(n_inner + n + 1):end, ik] = vec(GYdn[:, :, ik])
        end

        return nothing
    end

    return f, g!
end

"""
    disentangle(model; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Run disentangle on a `MagModel`.

Wannierize the up and down spin channels using the same `U` matrix.

# Arguments
- `model`: MagModel

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function disentangle(
    model::MagModel{T};
    λ::T=1.0,
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
    n_bands = model.up.n_bands
    n_wann = model.up.n_wann
    n_kpts = model.up.n_kpts

    @assert model.dn.n_bands == n_bands
    @assert model.dn.n_wann == n_wann
    @assert model.dn.n_kpts == n_kpts

    # TODO need QR orthogonalization rather than SVD to preserve the sparsity structure of Y?
    XYk_up_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (n_wann, n_wann), (n_bands, n_wann)
    )
    XYk_dn_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (n_wann, n_wann), (n_bands, n_wann)
    )
    n_inner = n_wann^2 + n_bands * n_wann
    XYkManif = Optim.ProductManifold(XYk_up_Manif, XYk_dn_Manif, (n_inner,), (n_inner,))
    XYManif = Optim.PowerManifold(XYkManif, (2 * n_inner,), (n_kpts,))

    Xup0, Yup0 = U_to_X_Y(model.up.U, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.U, model.dn.frozen_bands)
    # compact storage
    XYup0 = X_Y_to_XY(Xup0, Yup0)
    XYdn0 = X_Y_to_XY(Xdn0, Ydn0)
    XY0 = zeros(eltype(XYup0), (2 * n_inner, n_kpts))
    XY0[1:n_inner, :] = XYup0
    XY0[(n_inner + 1):end, :] = XYdn0

    # We have three storage formats:
    # (X, Y): n_wann * n_wann * n_kpts, n_bands * n_wann * n_kpts
    # U: n_bands * n_wann * n_kpts
    # XY: (n_wann * n_wann + n_bands * n_wann) * n_kpts
    f, g! = get_fg!_disentangle(model, n_bands, n_wann, n_kpts, λ)

    Ωup = omega(model.up)
    Ωdn = omega(model.dn)
    @info "Initial spread"
    show(Ωup)
    println("")
    show(Ωdn)
    println("\n")

    Oupdn = updn_overlap(model, model.up.U, model.dn.U)
    @info "Initial ↑↓ overlap"
    display(Oupdn)
    println("\n")

    Ωup = omega(model.up, Xup0, Yup0)
    Ωdn = omega(model.dn, Xdn0, Ydn0)
    @info "Initial spread (with states freezed)"
    show(Ωup)
    println("")
    show(Ωdn)
    println("\n")

    Oupdn = updn_overlap(model, X_Y_to_U(Xup0, Yup0), X_Y_to_U(Xdn0, Ydn0))
    @info "Initial ↑↓ overlap"
    display(Oupdn)
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
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, n_bands, n_wann)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, n_bands, n_wann)
    Uupmin = X_Y_to_U(Xupmin, Yupmin)
    Udnmin = X_Y_to_U(Xdnmin, Ydnmin)

    Ωup = omega(model.up, Uupmin)
    Ωdn = omega(model.dn, Udnmin)
    @info "Final spread"
    show(Ωup)
    println("")
    show(Ωdn)
    println("\n")

    Oupdn = updn_overlap(model, Uupmin, Udnmin)
    @info "Final ↑↓ overlap"
    display(Oupdn)
    println("\n")

    return Uupmin, Udnmin
end
