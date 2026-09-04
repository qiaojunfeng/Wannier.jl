module WannierManoptExt

using Wannier
using Wannier: Problem, Variance, ULayout, ManoptLBFGS, Model
using Wannier: LocalizationResult, LocalizationTraceEntry
using Wannier: n_bands, n_wannier, n_kpoints, _optimizer_callback
using Manopt
using Manifolds
using LinearAlgebra: norm

# -------------------------------------------------------------------------
# Variance + ULayout (isolated max_localize)
#
# Build a PowerManifold over Stiefel(n_bands, n_wannier; field = ℂ), wrap
# the fused (F, G, U) closure from `_optimizer_callback` as a Manopt cost and
# Riemannian gradient. Manopt does its own workspace for the quasi-Newton
# state; the Euclidean gradient buffer is preallocated here and shared
# between cost/grad calls to minimize allocations.
# -------------------------------------------------------------------------

function Wannier.solve(
        prob::Problem{<:Variance, <:Model, <:ULayout},
        solver::ManoptLBFGS;
        warmup::Bool = false,
    )
    model = prob.model
    nb = n_bands(model)
    nw = n_wannier(model)
    nk = n_kpoints(model)
    T = eltype(model.gauges)

    St = Stiefel(nb, nw, ℂ)
    M = PowerManifold(St, NestedPowerRepresentation(), nk)

    # Fused (F, G, U) closure shared with the OptimLBFGS path.
    fg! = _optimizer_callback(prob)

    # Reusable Euclidean buffers, copied in/out of the nested power repr.
    U3 = zeros(T, nb, nw, nk)
    G3 = zeros(T, nb, nw, nk)

    function _pack!(U3, p)
        @inbounds for ik in 1:nk
            view(U3, :, :, ik) .= p[ik]
        end
        return U3
    end

    cost(::AbstractManifold, p) = begin
        _pack!(U3, p)
        fg!(true, nothing, U3)::Float64
    end

    grad(Mp::AbstractManifold, p) = begin
        _pack!(U3, p)
        fg!(nothing, G3, U3)
        # Project Euclidean gradient onto tangent space of Stiefel at each k.
        gR = [project(St, p[ik], view(G3, :, :, ik)) for ik in 1:nk]
        return gR
    end

    # Initial point: current model gauges, copied into nested form.
    p0 = [collect(view(model.gauges, :, :, ik)) for ik in 1:nk]

    if warmup
        cost(M, p0)
        grad(M, p0)
    end

    stop = StopWhenGradientNormLess(solver.g_tol) |
        StopAfterIteration(solver.max_iter)

    record = solver.store_trace ? [:Iteration, :Cost, :GradientNorm] : missing
    debug = if solver.show_every > 0
        [
            :Iteration,
            " | f = ",
            :Cost,
            " | |g| = ",
            :GradientNorm,
            "\n",
            solver.show_every,
        ]
    else
        missing
    end
    start_time = time_ns()
    state = quasi_Newton(
        M, cost, grad, p0;
        memory_size = solver.memory_size,
        stopping_criterion = stop,
        record,
        debug,
        return_state = true,
    )
    elapsed_seconds = (time_ns() - start_time) / 1.0e9
    pmin = get_solver_result(state)

    Umin = zeros(T, nb, nw, nk)
    @inbounds for ik in 1:nk
        view(Umin, :, :, ik) .= pmin[ik]
    end
    trace = if solver.store_trace
        [
            LocalizationTraceEntry(Int(iteration), Float64(value), Float64(gradient_norm)) for
                (iteration, value, gradient_norm) in get_record(state)
        ]
    else
        LocalizationTraceEntry{Float64}[]
    end
    objective_value = Float64(cost(M, pmin))
    gradient_norm = Float64(norm(M, pmin, get_gradient(state)))
    converged = has_converged(state)
    termination_reason = converged ? :gradient_tolerance : :iteration_limit
    iterations = get_count(get_state(state), :Iterations)
    return LocalizationResult(
        Umin,
        objective_value,
        gradient_norm,
        iterations,
        converged,
        termination_reason,
        elapsed_seconds,
        trace,
    )
end

# Other Problem variants (XYLayout / ProductLayout / WLayout) fall through
# to the catch-all in src/localization/solver.jl, which emits a
# not-yet-implemented message. Add concrete methods here as they land.

end # module
