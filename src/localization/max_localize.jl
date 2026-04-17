using Optim: Optim

export localize_isolated_bands, max_localize

"""
    get_fg!_maxloc(model::Model)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_maxloc(terms::Tuple, model::Model)
    problem = LocalizationProblem(terms, model, :maxloc)
    return build_fg!(problem)
end

get_fg!_maxloc(model::Model) = get_fg!_maxloc((VarianceTerm(),), model)

"""
    max_localize(model; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Maximally localize spread functional w.r.t. all kpoints on a unitary matrix manifold.

# Arguments
- `model`: model

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function localize_isolated_bands(
    terms::Tuple,
        model::Model{T};
        f_tol::T = 1.0e-7,
        g_tol::T = 1.0e-5,
        max_iter::Int = 200,
        history_size::Int = 3,
    ) where {T <: Real}
    n_bands(model) != n_wannier(model) &&
        error("n_bands != n_wann, run instead disentanglement?")

    fg! = get_fg!_maxloc(terms, model)

    Ωⁱ = omega(model.kstencil, model.overlaps, model.gauges)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    Manif = manifold(UGauge(), model)

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    opt = Optim.optimize(
        Optim.only_fg!(fg!),
        model.gauges,
        meth(; manifold = Manif, linesearch = ls, m = history_size),
        Optim.Options(;
            show_trace = true,
            iterations = max_iter,
            f_tol = f_tol,
            g_tol = g_tol,
            allow_f_increases = true,
        ),
    )
    display(opt)

    Umin = Optim.minimizer(opt)

    Ωᶠ = omega(model.kstencil, model.overlaps, Umin)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Umin
end

localize_isolated_bands(model::Model; kwargs...) =
    localize_isolated_bands((VarianceTerm(),), model; kwargs...)

max_localize(terms::Tuple, model::Model{T}; kwargs...) where {T <: Real} =
    localize_isolated_bands(terms, model; kwargs...)

max_localize(model::Model; kwargs...) =
    solve!(Problem(Variance(), model), OptimLBFGS(; kwargs...))
