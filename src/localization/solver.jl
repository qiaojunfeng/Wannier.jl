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
            "The current extension covers Variance + ULayout (isolated " *
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

No rescaling happens here: what the objective produces is already the gradient
of the value this closure returns, in the convention Optim.jl consumes. See the
gradient-convention block at the top of `src/spread.jl`.
"""
function _optimizer_callback end

# One bridge for every (objective, model, layout): the objective computes in
# canonical coordinates, the layout converts. `WLayout` keeps a bespoke method
# below because its gradient is a contraction over kpoints rather than a repack.

_canonical_gradient(ws::Workspace) = ws.GU
_canonical_gradient(ws::SpinWorkspace) = (ws.up.GU, ws.dn.GU)

function _optimizer_callback(prob::Problem)
    obj, model, layout, ws = prob.objective, prob.model, prob.layout, prob.workspace
    GU = _canonical_gradient(ws)
    return function (F, G, x)
        U = assemble_gauge!(layout, x, model, ws)
        Ω = fg!(F, G === nothing ? nothing : GU, obj, U, model, ws)
        G === nothing || pullback_gradient!(G, layout, model, ws)
        return Ω
    end
end

# -------------------------------------------------------------------------
# solve! bindings
# -------------------------------------------------------------------------

function _run_optim(fg!, x0, man, solver::OptimLBFGS)
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

function solve!(prob::Problem, solver::OptimLBFGS)
    model, layout = prob.model, prob.layout
    opt = _run_optim(
        _optimizer_callback(prob), initial_parameters(layout, model), manifold(layout, model), solver
    )
    return finalize_result(layout, Optim.minimizer(opt), model)
end

solve!(prob::Problem; kwargs...) = solve!(prob, OptimLBFGS(; kwargs...))
