using Optim: Optim
using LinearAlgebra: I

"""
    AbstractLocalizationSolver

Solver-side contract. `solve!(prob::Problem, solver::AbstractLocalizationSolver)`
is the single entry point; solvers own the backend choice, tolerances,
iteration limits, linesearch, and history size.
"""
abstract type AbstractLocalizationSolver end

"""
    OptimLBFGS(; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3, linesearch=Optim.HagerZhang())

Optim.jl LBFGS driver.
"""
struct OptimLBFGS{LS} <: AbstractLocalizationSolver
    f_tol::Float64
    g_tol::Float64
    max_iter::Int
    history_size::Int
    linesearch::LS
end

function OptimLBFGS(;
        f_tol::Real = 1.0e-7,
        g_tol::Real = 1.0e-5,
        max_iter::Integer = 200,
        history_size::Integer = 3,
        linesearch = Optim.HagerZhang(),
    )
    return OptimLBFGS(Float64(f_tol), Float64(g_tol), Int(max_iter), Int(history_size), linesearch)
end

"""
    solve!(prob, solver)

Run `solver` against `prob`. Returns the optimized gauge (or gauge tuple
for `SpinModel`); does not mutate `prob.model`.
"""
function solve! end

# -------------------------------------------------------------------------
# Optim fused fg! closures per (Objective, Model, Layout)
# -------------------------------------------------------------------------

"""
Build an `Optim.only_fg!`-compatible `fg!(F, G, x)` closure for `prob`.
Writes gradients into `G` (the layout-native buffer Optim hands us) and
returns Ω when `F !== nothing`.
"""
function _make_optim_fg! end

# --- Variance + UGauge (isolated maxloc) ---

function _make_optim_fg!(prob::Problem{<:Variance, <:Model, <:UGauge})
    model = prob.model
    ws = prob.workspace
    kstencil = model.kstencil
    overlaps = model.overlaps
    return function fg!(F, G, U)
        compute_MU_UtMU!(ws, kstencil, overlaps, U)
        if G !== nothing
            omega_grad!(G, ws, kstencil, overlaps)
        end
        if F === nothing
            return nothing
        end
        return omega!(ws, kstencil, overlaps).Ω
    end
end

# --- Variance + XYGauge (disentangle) ---

function _make_optim_fg!(prob::Problem{<:Variance, <:Model, <:XYGauge})
    model = prob.model
    ws = prob.workspace
    kstencil = model.kstencil
    overlaps = model.overlaps
    frozen = model.frozen_bands
    nw = n_wannier(model)
    nk = n_kpoints(model)
    n = nw^2
    return function fg!(F, G, XY)
        X, Y = XY_to_X_Y!(ws.X, ws.Y, XY)
        U = X_Y_to_U!(ws.U, X, Y)
        compute_MU_UtMU!(ws, kstencil, overlaps, U)
        if G !== nothing
            omega_grad!(ws.G, ws, kstencil, overlaps)
            GX, GY = GU_to_GX_GY(ws.G, X, Y, frozen)
            @inbounds for ik in 1:nk
                gxk = view(GX, :, :, ik)
                gyk = view(GY, :, :, ik)
                for i in eachindex(gxk)
                    G[i, ik] = gxk[i]
                end
                for i in eachindex(gyk)
                    G[n + i, ik] = gyk[i]
                end
            end
        end
        if F === nothing
            return nothing
        end
        return omega!(ws, kstencil, overlaps).Ω
    end
end

# --- CenteredVariance + UGauge ---

function _make_optim_fg!(prob::Problem{<:CenteredVariance, <:Model, <:UGauge})
    obj = prob.objective
    model = prob.model
    ws = prob.workspace
    kstencil = model.kstencil
    overlaps = model.overlaps
    pen = center_penalty(obj.r0, obj.λ)
    return function fg!(F, G, U)
        compute_MU_UtMU!(ws, kstencil, overlaps, U)
        if G !== nothing
            omega_grad!(pen, G, ws, kstencil, overlaps)
        end
        Ωbase = omega!(ws, kstencil, overlaps)
        if F === nothing
            return nothing
        end
        return omega_center(Ωbase; r₀ = obj.r0, λ = obj.λ).Ωt
    end
end

# --- CenteredVariance + XYGauge ---

function _make_optim_fg!(prob::Problem{<:CenteredVariance, <:Model, <:XYGauge})
    obj = prob.objective
    model = prob.model
    ws = prob.workspace
    kstencil = model.kstencil
    overlaps = model.overlaps
    frozen = model.frozen_bands
    pen = center_penalty(obj.r0, obj.λ)
    nw = n_wannier(model)
    nk = n_kpoints(model)
    n = nw^2
    return function fg!(F, G, XY)
        X, Y = XY_to_X_Y!(ws.X, ws.Y, XY)
        U = X_Y_to_U!(ws.U, X, Y)
        compute_MU_UtMU!(ws, kstencil, overlaps, U)
        if G !== nothing
            omega_grad!(pen, ws.G, ws, kstencil, overlaps)
            GX, GY = GU_to_GX_GY(ws.G, X, Y, frozen)
            @inbounds for ik in 1:nk
                gxk = view(GX, :, :, ik)
                gyk = view(GY, :, :, ik)
                for i in eachindex(gxk)
                    G[i, ik] = gxk[i]
                end
                for i in eachindex(gyk)
                    G[n + i, ik] = gyk[i]
                end
            end
        end
        Ωbase = omega!(ws, kstencil, overlaps)
        if F === nothing
            return nothing
        end
        return omega_center(Ωbase; r₀ = obj.r0, λ = obj.λ).Ωt
    end
end

# --- CoOptVariance + SpinModel + ProductLayout ---

function _make_optim_fg!(prob::Problem{<:CoOptVariance, <:SpinModel, <:ProductLayout})
    model = prob.model
    ws = prob.workspace
    λ = prob.objective.λs
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)
    n_inner = nw^2 + nb * nw
    n = nw^2
    return function fg!(F, G, XY)
        XYr = reshape(XY, (2 * n_inner, nk))
        XYup = @view XYr[1:n_inner, :]
        XYdn = @view XYr[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)

        Ωtot = nothing
        if F !== nothing
            Ωup = omega(model.up.kstencil, model.up.overlaps, Xup, Yup).Ω
            Ωdn = omega(model.dn.kstencil, model.dn.overlaps, Xdn, Ydn).Ω
            Ωupdn = λ == 0 ? 0.0 :
                omega_updn(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
            Ωtot = Ωup + Ωdn + λ * Ωupdn
        end

        if G !== nothing
            GXup, GYup = omega_grad(
                model.up.kstencil, model.up.overlaps, Xup, Yup, model.up.frozen_bands
            )
            GXdn, GYdn = omega_grad(
                model.dn.kstencil, model.dn.overlaps, Xdn, Ydn, model.dn.frozen_bands
            )
            if λ != 0
                GOXup, GOYup, GOXdn, GOYdn = omega_updn_grad(model, Xup, Yup, Xdn, Ydn)
                GXup += λ * GOXup
                GYup += λ * GOYup
                GXdn += λ * GOXdn
                GYdn += λ * GOYdn
            end
            for ik in 1:nk
                for (i, v) in enumerate(view(GXup, :, :, ik)); G[i, ik] = v; end
                for (i, v) in enumerate(view(GYup, :, :, ik)); G[n + i, ik] = v; end
                for (i, v) in enumerate(view(GXdn, :, :, ik)); G[n_inner + i, ik] = v; end
                for (i, v) in enumerate(view(GYdn, :, :, ik)); G[n_inner + n + i, ik] = v; end
            end
        end

        return F === nothing ? nothing : Ωtot
    end
end

# --- CenteredCoOptVariance + SpinModel + ProductLayout ---

function _make_optim_fg!(prob::Problem{<:CenteredCoOptVariance, <:SpinModel, <:ProductLayout})
    model = prob.model
    obj = prob.objective
    λs = obj.λs
    pen = center_penalty(obj.r0, obj.λ)
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)
    n_inner = nw^2 + nb * nw
    n = nw^2

    function f(XY)
        XYr = reshape(XY, (2 * n_inner, nk))
        XYup = @view XYr[1:n_inner, :]
        XYdn = @view XYr[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        Ωup = omega_center(
            omega(model.up.kstencil, model.up.overlaps, Xup, Yup); r₀ = obj.r0, λ = obj.λ
        ).Ωt
        Ωdn = omega_center(
            omega(model.dn.kstencil, model.dn.overlaps, Xdn, Ydn); r₀ = obj.r0, λ = obj.λ
        ).Ωt
        Ωupdn = λs == 0 ? 0.0 :
            omega_updn(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
        return Ωup + Ωdn + λs * Ωupdn
    end

    function g!(G, XY)
        XYr = reshape(XY, (2 * n_inner, nk))
        XYup = @view XYr[1:n_inner, :]
        XYdn = @view XYr[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        GXup, GYup = omega_grad(
            pen, model.up.kstencil, model.up.overlaps, Xup, Yup, model.up.frozen_bands
        )
        GXdn, GYdn = omega_grad(
            pen, model.dn.kstencil, model.dn.overlaps, Xdn, Ydn, model.dn.frozen_bands
        )
        if λs != 0
            GOXup, GOYup, GOXdn, GOYdn = omega_updn_grad(model, Xup, Yup, Xdn, Ydn)
            GXup += λs * GOXup
            GYup += λs * GOYup
            GXdn += λs * GOXdn
            GYdn += λs * GOYdn
        end
        for ik in 1:nk
            G[1:n, ik] = vec(view(GXup, :, :, ik))
            G[(n + 1):n_inner, ik] = vec(view(GYup, :, :, ik))
            G[(n_inner + 1):(n_inner + n), ik] = vec(view(GXdn, :, :, ik))
            G[(n_inner + n + 1):end, ik] = vec(view(GYdn, :, :, ik))
        end
        return nothing
    end

    return f, g!
end

# -------------------------------------------------------------------------
# solve! bindings
# -------------------------------------------------------------------------

function _run_optim_fg!(fg!, x0, man, solver::OptimLBFGS)
    return Optim.optimize(
        Optim.only_fg!(fg!),
        x0,
        Optim.LBFGS(; manifold = man, linesearch = solver.linesearch, m = solver.history_size),
        Optim.Options(;
            f_tol = solver.f_tol,
            g_tol = solver.g_tol,
            iterations = solver.max_iter,
            allow_f_increases = true,
            show_trace = true,
        ),
    )
end

function _run_optim_fg!(f::Function, g!::Function, x0, man, solver::OptimLBFGS)
    return Optim.optimize(
        f, g!, x0,
        Optim.LBFGS(; manifold = man, linesearch = solver.linesearch, m = solver.history_size),
        Optim.Options(;
            f_tol = solver.f_tol,
            g_tol = solver.g_tol,
            iterations = solver.max_iter,
            allow_f_increases = true,
            show_trace = true,
        ),
    )
end

# Variance on Model (UGauge or XYGauge)
function solve!(prob::Problem{<:Variance, <:Model, <:Union{UGauge, XYGauge}}, solver::OptimLBFGS)
    model = prob.model
    layout = prob.layout
    fg! = _make_optim_fg!(prob)
    x0 = if layout isa UGauge
        copy(model.gauges)
    else
        X0, Y0 = U_to_X_Y(model.gauges, model.frozen_bands)
        X_Y_to_XY(X0, Y0)
    end
    man = manifold(layout, model)
    opt = _run_optim_fg!(fg!, x0, man, solver)
    xmin = Optim.minimizer(opt)
    return layout isa UGauge ? xmin :
        let (Xm, Ym) = XY_to_X_Y(xmin, n_bands(model), n_wannier(model))
            X_Y_to_U(Xm, Ym)
        end
end

# CenteredVariance on Model
function solve!(prob::Problem{<:CenteredVariance, <:Model}, solver::OptimLBFGS)
    model = prob.model
    layout = prob.layout
    fg! = _make_optim_fg!(prob)
    x0 = if layout isa UGauge
        copy(model.gauges)
    elseif layout isa XYGauge
        X0, Y0 = U_to_X_Y(model.gauges, model.frozen_bands)
        X_Y_to_XY(X0, Y0)
    else
        error("solve!(CenteredVariance, ::$(typeof(layout))) not supported yet")
    end
    man = manifold(layout, model)
    opt = _run_optim_fg!(fg!, x0, man, solver)
    xmin = Optim.minimizer(opt)
    return layout isa UGauge ? xmin :
        let (Xm, Ym) = XY_to_X_Y(xmin, n_bands(model), n_wannier(model))
            X_Y_to_U(Xm, Ym)
        end
end

# Variance + WLayout (opt_rotate) — keeps the rotation-specific fg!
function solve!(prob::Problem{<:Variance, <:Model, <:WLayout}, solver::OptimLBFGS)
    model = prob.model
    n_wann = n_wannier(model)
    n_bands(model) == n_wann || error("n_bands != n_wann, run disentanglement instead")
    model2 = deepcopy(model)
    model2.overlaps .= transform_gauge(model2.overlaps, model2.kstencil.kpb_k, model2.gauges)
    model2.gauges  .= identity_gauge(eltype(model2.gauges), n_kpoints(model2), n_wann)
    f, g! = get_fg!_rotate(model2)
    man   = manifold(WLayout(), model2)
    W0    = Matrix{eltype(model2.gauges)}(I, n_wann, n_wann)
    opt   = _run_optim_fg!(f, g!, W0, man, solver)
    return Optim.minimizer(opt)
end

# CoOptVariance on SpinModel
function solve!(prob::Problem{<:CoOptVariance, <:SpinModel, <:ProductLayout}, solver::OptimLBFGS)
    model = prob.model
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    n_inner = nw^2 + nb * nw
    fg! = _make_optim_fg!(prob)
    Xup0, Yup0 = U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
    XY0 = vcat(X_Y_to_XY(Xup0, Yup0), X_Y_to_XY(Xdn0, Ydn0))
    man = manifold(prob.layout, model)
    opt = _run_optim_fg!(fg!, XY0, man, solver)
    XYmin = Optim.minimizer(opt)
    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nb, nw)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nb, nw)
    return X_Y_to_U(Xupmin, Yupmin), X_Y_to_U(Xdnmin, Ydnmin)
end

# CenteredCoOptVariance on SpinModel
function solve!(prob::Problem{<:CenteredCoOptVariance, <:SpinModel, <:ProductLayout}, solver::OptimLBFGS)
    model = prob.model
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    n_inner = nw^2 + nb * nw
    f, g! = _make_optim_fg!(prob)
    Xup0, Yup0 = U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
    XY0 = vcat(X_Y_to_XY(Xup0, Yup0), X_Y_to_XY(Xdn0, Ydn0))
    man = manifold(prob.layout, model)
    opt = _run_optim_fg!(f, g!, XY0, man, solver)
    XYmin = Optim.minimizer(opt)
    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nb, nw)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nb, nw)
    return X_Y_to_U(Xupmin, Yupmin), X_Y_to_U(Xdnmin, Ydnmin)
end

solve!(prob::Problem; kwargs...) = solve!(prob, OptimLBFGS(; kwargs...))
