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

Optim.jl LBFGS driver. Mirrors the tolerances / options of the legacy
`disentangle` / `max_localize` entry points so migration is a pure
call-site rewrite.
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

# Concrete bindings for `Variance` on a `Model` live here. Other Objective
# subtypes (CenteredVariance, CoOptVariance, CenteredCoOptVariance) migrate
# during Phase 5 (commits T–W); until then the legacy entry points keep
# driving Optim directly.

function solve!(prob::Problem{<:Variance, <:Model}, solver::OptimLBFGS)
    model = prob.model
    layout = prob.layout
    if layout isa UGauge
        legacy = LocalizationProblem((VarianceTerm(),), model, :maxloc)
        fg! = _build_fg_maxloc(legacy)
        x0 = copy(model.gauges)
    elseif layout isa XYGauge
        legacy = LocalizationProblem((VarianceTerm(),), model, :disentangle)
        fg! = _build_fg_disentangle(legacy)
        X0, Y0 = U_to_X_Y(model.gauges, model.frozen_bands)
        x0 = X_Y_to_XY(X0, Y0)
    else
        error("solve!(Variance, ::$(typeof(layout))) not supported yet")
    end

    man = manifold(layout, model)
    opt = Optim.optimize(
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
    xmin = Optim.minimizer(opt)

    if layout isa UGauge
        return xmin
    else
        Xmin, Ymin = XY_to_X_Y(xmin, n_bands(model), n_wannier(model))
        return X_Y_to_U(Xmin, Ymin)
    end
end

function solve!(prob::Problem{<:Variance, <:Model, <:WLayout}, solver::OptimLBFGS)
    model = prob.model
    n_wann = n_wannier(model)
    n_bands(model) == n_wann || error("n_bands != n_wann, run disentanglement instead")

    # Work on a copy with transformed overlaps + identity gauges so the
    # optimization variable is purely the W rotation.
    model2 = deepcopy(model)
    model2.overlaps .= transform_gauge(model2.overlaps, model2.kstencil.kpb_k, model2.gauges)
    model2.gauges  .= identity_gauge(eltype(model2.gauges), n_kpoints(model2), n_wann)

    f, g! = get_fg!_rotate(model2)
    man   = manifold(WLayout(), model2)
    W0    = Matrix{eltype(model2.gauges)}(I, n_wann, n_wann)

    opt = Optim.optimize(
        f, g!, W0,
        Optim.LBFGS(; manifold = man, linesearch = solver.linesearch, m = solver.history_size),
        Optim.Options(;
            f_tol = solver.f_tol,
            g_tol = solver.g_tol,
            iterations = solver.max_iter,
            allow_f_increases = true,
            show_trace = true,
        ),
    )
    return Optim.minimizer(opt)
end

function solve!(
        prob::Problem{<:CoOptVariance, <:SpinModel, <:ProductLayout}, solver::OptimLBFGS
    )
    model = prob.model
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    n_inner = nw^2 + nb * nw

    legacy = LocalizationProblem(
        (VarianceTerm(),), model, :mag_disentangle; lambda = prob.objective.λs
    )
    _, _, fg! = _build_fg_mag_disentangle(legacy)

    Xup0, Yup0 = U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
    XY0 = vcat(X_Y_to_XY(Xup0, Yup0), X_Y_to_XY(Xdn0, Ydn0))

    man = manifold(prob.layout, model)
    opt = Optim.optimize(
        Optim.only_fg!(fg!),
        XY0,
        Optim.LBFGS(; manifold = man, linesearch = solver.linesearch, m = solver.history_size),
        Optim.Options(;
            f_tol = solver.f_tol,
            g_tol = solver.g_tol,
            iterations = solver.max_iter,
            allow_f_increases = true,
            show_trace = true,
        ),
    )
    XYmin = Optim.minimizer(opt)
    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nb, nw)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nb, nw)
    return X_Y_to_U(Xupmin, Yupmin), X_Y_to_U(Xdnmin, Ydnmin)
end

function solve!(
        prob::Problem{<:CenteredCoOptVariance, <:SpinModel, <:ProductLayout}, solver::OptimLBFGS
    )
    model = prob.model
    obj = prob.objective
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    n_inner = nw^2 + nb * nw

    terms = (VarianceTerm(), CenterConstraintTerm(obj.r0, obj.λ))
    legacy = LocalizationProblem(terms, model, :mag_disentangle_center; lambda = obj.λs)
    f, g! = _build_fg_mag_disentangle_center(legacy)

    Xup0, Yup0 = U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
    XY0 = vcat(X_Y_to_XY(Xup0, Yup0), X_Y_to_XY(Xdn0, Ydn0))

    man = manifold(prob.layout, model)
    opt = Optim.optimize(
        f, g!, XY0,
        Optim.LBFGS(; manifold = man, linesearch = solver.linesearch, m = solver.history_size),
        Optim.Options(;
            f_tol = solver.f_tol,
            g_tol = solver.g_tol,
            iterations = solver.max_iter,
            allow_f_increases = true,
            show_trace = true,
        ),
    )
    XYmin = Optim.minimizer(opt)
    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nb, nw)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nb, nw)
    return X_Y_to_U(Xupmin, Yupmin), X_Y_to_U(Xdnmin, Ydnmin)
end

function solve!(prob::Problem{<:CenteredVariance, <:Model}, solver::OptimLBFGS)
    model = prob.model
    obj = prob.objective
    layout = prob.layout
    terms = (VarianceTerm(), CenterConstraintTerm(obj.r0, obj.λ))

    if layout isa UGauge
        legacy = LocalizationProblem(terms, model, :maxloc)
        fg! = _build_fg_maxloc(legacy)
        x0 = copy(model.gauges)
    elseif layout isa XYGauge
        legacy = LocalizationProblem(terms, model, :disentangle)
        fg! = _build_fg_disentangle(legacy)
        X0, Y0 = U_to_X_Y(model.gauges, model.frozen_bands)
        x0 = X_Y_to_XY(X0, Y0)
    else
        error("solve!(CenteredVariance, ::$(typeof(layout))) not supported yet")
    end

    man = manifold(layout, model)
    opt = Optim.optimize(
        Optim.only_fg!(fg!), x0,
        Optim.LBFGS(; manifold = man, linesearch = solver.linesearch, m = solver.history_size),
        Optim.Options(;
            f_tol = solver.f_tol,
            g_tol = solver.g_tol,
            iterations = solver.max_iter,
            allow_f_increases = true,
            show_trace = true,
        ),
    )
    xmin = Optim.minimizer(opt)

    if layout isa UGauge
        return xmin
    else
        Xmin, Ymin = XY_to_X_Y(xmin, n_bands(model), n_wannier(model))
        return X_Y_to_U(Xmin, Ymin)
    end
end

"""
    solve!(prob; kwargs...) = solve!(prob, OptimLBFGS(; kwargs...))

Convenience wrapper for callers that don't care about solver backend.
"""
solve!(prob::Problem; kwargs...) = solve!(prob, OptimLBFGS(; kwargs...))
