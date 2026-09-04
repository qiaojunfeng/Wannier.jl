using Optim: Optim
using LinearAlgebra: I

export solve, LocalizationResult, LocalizationTraceEntry, OptimLBFGS, ManoptLBFGS

"""
    LocalizationTraceEntry

Backend-independent diagnostics for one localization iteration.
"""
struct LocalizationTraceEntry{T <: Real}
    iteration::Int
    objective_value::T
    gradient_norm::T
end

"""
    LocalizationResult

Result returned by [`solve`](@ref). `solution` has the same form that
[`localize`](@ref) returns: a gauge, a gauge tuple for a spin model, or the
`(U_fbz, U_ibz)` pair for a symmetry-constrained model.

`termination_reason` is one of `:gradient_tolerance`, `:objective_tolerance`,
`:parameter_tolerance`, `:iteration_limit`, `:line_search_failure`, or
`:unknown`. The optional `trace` contains only scalar, backend-independent
diagnostics; backend-specific optimizer state is deliberately not retained.
"""
struct LocalizationResult{S, T <: Real}
    solution::S
    objective_value::T
    gradient_norm::T
    iterations::Int
    converged::Bool
    termination_reason::Symbol
    elapsed_seconds::Float64
    trace::Vector{LocalizationTraceEntry{T}}
end

function Base.show(io::IO, ::MIME"text/plain", result::LocalizationResult)
    println(io, "LocalizationResult:")
    println(io, "  converged           = ", result.converged)
    println(io, "  termination_reason  = ", result.termination_reason)
    println(io, "  objective_value     = ", result.objective_value)
    println(io, "  gradient_norm       = ", result.gradient_norm)
    println(io, "  iterations          = ", result.iterations)
    println(io, "  elapsed_seconds     = ", result.elapsed_seconds)
    return print(io, "  trace entries       = ", length(result.trace))
end

"""
    AbstractLocalizationSolver

Solver-side contract. `solve(prob::Problem, solver::AbstractLocalizationSolver)`
is the single entry point; solvers own the backend choice, tolerances,
iteration limits, linesearch, and history size.
"""
abstract type AbstractLocalizationSolver end

"""
    OptimLBFGS(; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3,
                linesearch=Optim.HagerZhang(), store_trace=false, show_every=0)

Optim.jl LBFGS driver.
"""
struct OptimLBFGS{LS} <: AbstractLocalizationSolver
    f_tol::Float64
    g_tol::Float64
    max_iter::Int
    history_size::Int
    linesearch::LS
    store_trace::Bool
    show_every::Int
end

function OptimLBFGS(;
        f_tol::Real = 1.0e-7,
        g_tol::Real = 1.0e-5,
        max_iter::Integer = 200,
        history_size::Integer = 3,
        linesearch = Optim.HagerZhang(),
        store_trace::Bool = false,
        show_every::Integer = 0,
    )
    show_every >= 0 || throw(ArgumentError("show_every must be nonnegative"))
    return OptimLBFGS(
        Float64(f_tol),
        Float64(g_tol),
        Int(max_iter),
        Int(history_size),
        linesearch,
        store_trace,
        Int(show_every),
    )
end

"""
    ManoptLBFGS(; g_tol=1e-5, max_iter=200, memory_size=3,
                 store_trace=false, show_every=0)

Manopt.jl L-BFGS driver on the Stiefel manifold. Requires loading
`Manopt` and `Manifolds` — the actual `solve` method lives in the
`WannierManoptExt` package extension and is only compiled once those
packages are `using`-ed.

# Keyword arguments
- `g_tol`: Riemannian gradient-norm stopping tolerance.
- `max_iter`: iteration cap.
- `memory_size`: L-BFGS history length.
- `store_trace`: retain scalar per-iteration diagnostics in the result.
- `show_every`: print every `n`th iteration, or remain quiet when zero.

See also [`OptimLBFGS`](@ref) for the default Optim.jl backend.
"""
struct ManoptLBFGS <: AbstractLocalizationSolver
    g_tol::Float64
    max_iter::Int
    memory_size::Int
    store_trace::Bool
    show_every::Int
end

function ManoptLBFGS(;
        g_tol::Real = 1.0e-5,
        max_iter::Integer = 200,
        memory_size::Integer = 3,
        store_trace::Bool = false,
        show_every::Integer = 0,
    )
    show_every >= 0 || throw(ArgumentError("show_every must be nonnegative"))
    return ManoptLBFGS(
        Float64(g_tol), Int(max_iter), Int(memory_size), store_trace, Int(show_every)
    )
end

# Fallback. The extension in ext/WannierManoptExt registers concrete
# methods on supported (Objective, Model, Layout) combinations; anything
# not overridden lands here. Works whether the extension is loaded or not.
function solve(::Problem, ::ManoptLBFGS; warmup::Bool = false)
    return error(
        "ManoptLBFGS is not available for this (Objective, Model, Layout). " *
            "The current extension covers Variance + ULayout (isolated " *
            "max_localize); other variants still route through OptimLBFGS. " *
            "If Manopt / Manifolds are not yet loaded, `using Manopt, " *
            "Manifolds` will activate the extension."
    )
end

"""
    solve(prob, solver; warmup=false)

Run `solver` against `prob` and return a [`LocalizationResult`](@ref). This
does not mutate `prob.model`. With `warmup=true`, perform one value-and-gradient
evaluation before starting the solver timer; this is useful when reporting
timings from a freshly started Julia process.
"""
function solve end

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
# solve bindings
# -------------------------------------------------------------------------

function _run_optim(fg!, x0, man, solver::OptimLBFGS)
    return Optim.optimize(
        Optim.only_fg!(fg!),
        x0,
        Optim.LBFGS(; manifold = man, linesearch = solver.linesearch, m = solver.history_size),
        Optim.Options(;
            f_reltol = solver.f_tol,
            g_tol = solver.g_tol,
            iterations = solver.max_iter,
            allow_f_increases = true,
            store_trace = solver.store_trace,
            show_trace = solver.show_every > 0,
            show_every = max(1, solver.show_every),
        ),
    )
end

function _optim_termination_reason(opt)
    Optim.g_converged(opt) && return :gradient_tolerance
    Optim.f_converged(opt) && return :objective_tolerance
    Optim.x_converged(opt) && return :parameter_tolerance
    Optim.iteration_limit_reached(opt) && return :iteration_limit
    Optim.termination_code(opt) == Optim.TerminationCode.FailedLinesearch &&
        return :line_search_failure
    return :unknown
end

function _optim_trace(opt, store_trace)
    store_trace || return LocalizationTraceEntry{Float64}[]
    return [
        LocalizationTraceEntry(
                state.iteration, Float64(state.value), Float64(state.g_norm)
            ) for state in Optim.trace(opt)
    ]
end

function solve(prob::Problem, solver::OptimLBFGS; warmup::Bool = false)
    model, layout = prob.model, prob.layout
    callback = _optimizer_callback(prob)
    initial = initial_parameters(layout, model)
    if warmup
        callback(true, similar(initial), initial)
    end
    start_time = time_ns()
    opt = _run_optim(callback, initial, manifold(layout, model), solver)
    elapsed_seconds = (time_ns() - start_time) / 1.0e9
    solution = finalize_result(layout, Optim.minimizer(opt), model)
    return LocalizationResult(
        solution,
        Float64(Optim.minimum(opt)),
        Float64(Optim.g_residual(opt)),
        Optim.iterations(opt),
        Optim.converged(opt),
        _optim_termination_reason(opt),
        elapsed_seconds,
        _optim_trace(opt, solver.store_trace),
    )
end

function solve(prob::Problem; warmup::Bool = false, kwargs...)
    return solve(prob, OptimLBFGS(; kwargs...); warmup)
end
