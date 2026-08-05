using Optim: Optim
using LinearAlgebra: I

export solve!, OptimLBFGS, ManoptLBFGS

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
    ManoptLBFGS(; g_tol=1e-5, max_iter=200, memory_size=3)

Manopt.jl L-BFGS driver on the Stiefel manifold. Requires loading
`Manopt` and `Manifolds` — the actual `solve!` method lives in the
`WannierManoptExt` package extension and is only compiled once those
packages are `using`-ed.

# Keyword arguments
- `g_tol`: Riemannian gradient-norm stopping tolerance.
- `max_iter`: iteration cap.
- `memory_size`: L-BFGS history length.

See also [`OptimLBFGS`](@ref) for the default Optim.jl backend.
"""
struct ManoptLBFGS <: AbstractLocalizationSolver
    g_tol::Float64
    max_iter::Int
    memory_size::Int
end

function ManoptLBFGS(; g_tol::Real = 1.0e-5, max_iter::Integer = 200, memory_size::Integer = 3)
    return ManoptLBFGS(Float64(g_tol), Int(max_iter), Int(memory_size))
end

# Fallback. The extension in ext/WannierManoptExt registers concrete
# methods on supported (Objective, Model, Layout) combinations; anything
# not overridden lands here. Works whether the extension is loaded or not.
function solve!(::Problem, ::ManoptLBFGS)
    return error(
        "ManoptLBFGS is not available for this (Objective, Model, Layout). " *
            "The current extension covers Variance + UGauge (isolated " *
            "max_localize); other variants still route through OptimLBFGS. " *
            "If Manopt / Manifolds are not yet loaded, `using Manopt, " *
            "Manifolds` will activate the extension."
    )
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

## Gradient convention

Internally, `omega_grad!` and related kernel functions compute gradients under
the physics convention: `df(x) = 2 Re⟨∇f, dx⟩`. Optim.jl expects the other
convention: `df = Re⟨∇f, dx⟩`. A factor of ½ should be applied to gradients
handed to Optim at this boundary. Currently the line-search has absorbed this
factor; an explicit rescale is deferred to improve robustness.

See also [`spread.jl`](@ref) gradient-convention block for full context.
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
    ws = prob.workspace::SpinWorkspace
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
        Xup, Yup = XY_to_X_Y!(ws.up.X, ws.up.Y, XYup)
        Xdn, Ydn = XY_to_X_Y!(ws.dn.X, ws.dn.Y, XYdn)
        Uup = X_Y_to_U!(ws.up.U, Xup, Yup)
        Udn = X_Y_to_U!(ws.dn.U, Xdn, Ydn)

        compute_MU_UtMU!(ws.up, model.up.kstencil, model.up.overlaps, Uup)
        compute_MU_UtMU!(ws.dn, model.dn.kstencil, model.dn.overlaps, Udn)

        if G !== nothing
            omega_grad!(ws.up.G, ws.up, model.up.kstencil, model.up.overlaps)
            omega_grad!(ws.dn.G, ws.dn, model.dn.kstencil, model.dn.overlaps)
            GXup, GYup = GU_to_GX_GY(ws.up.G, Xup, Yup, model.up.frozen_bands)
            GXdn, GYdn = GU_to_GX_GY(ws.dn.G, Xdn, Ydn, model.dn.frozen_bands)
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

        if F === nothing
            return nothing
        end
        Ωup = omega!(ws.up, model.up.kstencil, model.up.overlaps).Ω
        Ωdn = omega!(ws.dn, model.dn.kstencil, model.dn.overlaps).Ω
        Ωupdn = λ == 0 ? 0.0 : omega_updn(model, Uup, Udn)
        return Ωup + Ωdn + λ * Ωupdn
    end
end

# --- CenteredCoOptVariance + SpinModel + ProductLayout ---

function _make_optim_fg!(prob::Problem{<:CenteredCoOptVariance, <:SpinModel, <:ProductLayout})
    model = prob.model
    obj = prob.objective
    ws = prob.workspace::SpinWorkspace
    λs = obj.λs
    pen = center_penalty(obj.r0, obj.λ)
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)
    n_inner = nw^2 + nb * nw
    n = nw^2

    function _decode!(XY)
        XYr = reshape(XY, (2 * n_inner, nk))
        Xup, Yup = XY_to_X_Y!(ws.up.X, ws.up.Y, @view XYr[1:n_inner, :])
        Xdn, Ydn = XY_to_X_Y!(ws.dn.X, ws.dn.Y, @view XYr[(n_inner + 1):end, :])
        Uup = X_Y_to_U!(ws.up.U, Xup, Yup)
        Udn = X_Y_to_U!(ws.dn.U, Xdn, Ydn)
        compute_MU_UtMU!(ws.up, model.up.kstencil, model.up.overlaps, Uup)
        compute_MU_UtMU!(ws.dn, model.dn.kstencil, model.dn.overlaps, Udn)
        return Xup, Yup, Xdn, Ydn, Uup, Udn
    end

    function f(XY)
        _, _, _, _, Uup, Udn = _decode!(XY)
        Ωup = omega_center(
            omega!(ws.up, model.up.kstencil, model.up.overlaps); r₀ = obj.r0, λ = obj.λ
        ).Ωt
        Ωdn = omega_center(
            omega!(ws.dn, model.dn.kstencil, model.dn.overlaps); r₀ = obj.r0, λ = obj.λ
        ).Ωt
        Ωupdn = λs == 0 ? 0.0 : omega_updn(model, Uup, Udn)
        return Ωup + Ωdn + λs * Ωupdn
    end

    function g!(G, XY)
        Xup, Yup, Xdn, Ydn, _, _ = _decode!(XY)
        omega_grad!(pen, ws.up.G, ws.up, model.up.kstencil, model.up.overlaps)
        omega_grad!(pen, ws.dn.G, ws.dn, model.dn.kstencil, model.dn.overlaps)
        GXup, GYup = GU_to_GX_GY(ws.up.G, Xup, Yup, model.up.frozen_bands)
        GXdn, GYdn = GU_to_GX_GY(ws.dn.G, Xdn, Ydn, model.dn.frozen_bands)
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

# --- Variance + WLayout (opt_rotate) — single-W rotation manifold ---
#
# The optimization variable is the single rotation matrix `W`; `prob.model`
# is assumed to be in W-rotation-ready form (identity gauges, transformed
# overlaps). `solve!` below performs that transform once before constructing
# the closure.

function _make_optim_fg!(prob::Problem{<:Variance, <:Model, <:WLayout})
    model = prob.model
    M = model.overlaps
    bvectors = model.kstencil
    kpb_k = bvectors.kpb_k
    kpb_G = bvectors.kpb_G
    wb = bvectors.bweights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints
    n_bvecs = size(M, 3)
    n_kpts = n_kpoints(model)

    return function fg!(F, G, W)
        n_wann = size(W, 1)
        UW = merge_gauge(model.gauges, W)

        if G !== nothing
            T = zeros(eltype(W), n_wann, n_wann)
            MWᵏᵇ = zeros(eltype(W), n_wann, n_wann)
            Nᵏᵇ = zeros(eltype(W), n_wann, n_wann)
            r = center(bvectors, M, UW)

            fill!(G, 0.0)
            for ik in 1:n_kpts
                for ib in 1:n_bvecs
                    ikpb = kpb_k[ib, ik]
                    MWᵏᵇ .= view(M, :, :, ib, ik) * W
                    Nᵏᵇ .= W' * MWᵏᵇ
                    b = recip_lattice * (kpoints[ikpb] + kpb_G[ib, ik] - kpoints[ik])

                    q = imaglog.(diag(Nᵏᵇ))
                    for iw in 1:n_wann
                        q[iw] += r[iw] ⋅ b
                    end

                    for n in 1:n_wann
                        abs(Nᵏᵇ[n, n]) < 1.0e-10 && error("Nᵏᵇ too small! $ik -> $ikpb, $Nᵏᵇ")
                        t = -conj(Nᵏᵇ[n, n]) - im * q[n] / Nᵏᵇ[n, n]
                        for m in 1:n_wann
                            T[m, n] = t * MWᵏᵇ[m, n]
                        end
                    end

                    G .+= 4 * wb[ib] * T
                end
            end
            G ./= n_kpts
        end

        return F === nothing ? nothing : spread(model.kstencil, model.overlaps, UW).Ω
    end
end

function solve!(prob::Problem{<:Variance, <:Model, <:WLayout}, solver::OptimLBFGS)
    model = prob.model
    n_wann = n_wannier(model)
    n_bands(model) == n_wann || error("n_bands != n_wann, run disentanglement instead")
    model2 = deepcopy(model)
    model2.overlaps .= transform_gauge(model2.overlaps, model2.kstencil.kpb_k, model2.gauges)
    model2.gauges  .= identity_gauge(eltype(model2.gauges), n_kpoints(model2), n_wann)
    fg! = _make_optim_fg!(Problem(prob.objective, model2, prob.layout))
    man = manifold(WLayout(), model2)
    W0 = Matrix{eltype(model2.gauges)}(I, n_wann, n_wann)
    opt = _run_optim_fg!(fg!, W0, man, solver)
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
