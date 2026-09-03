module WannierManoptExt

using Wannier
using Wannier: Problem, Variance, ULayout, ManoptLBFGS, Model
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

function Wannier.solve!(
        prob::Problem{<:Variance, <:Model, <:ULayout}, solver::ManoptLBFGS
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

    stop = StopWhenGradientNormLess(solver.g_tol) |
        StopAfterIteration(solver.max_iter)

    pmin = quasi_Newton(
        M, cost, grad, p0;
        memory_size = solver.memory_size,
        stopping_criterion = stop,
    )

    Umin = zeros(T, nb, nw, nk)
    @inbounds for ik in 1:nk
        view(Umin, :, :, ik) .= pmin[ik]
    end
    return Umin
end

# Other Problem variants (XYLayout / ProductLayout / WLayout) fall through
# to the catch-all in src/localization/solver.jl, which emits a
# not-yet-implemented message. Add concrete methods here as they land.

end # module
