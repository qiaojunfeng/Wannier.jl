using Optim: Optim

export localize_isolated_bands, max_localize

"""
    get_fg!_maxloc(model::Model)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_maxloc(p::AbstractPenalty, model::Model)
    problem = LocalizationProblem(p, model, :maxloc)
    return build_fg!(problem)
end

get_fg!_maxloc(model::Model) = get_fg!_maxloc(SpreadPenalty(), model)

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
        p::AbstractPenalty,
        model::Model{T};
        f_tol::T = 1.0e-7,
        g_tol::T = 1.0e-5,
        max_iter::Int = 200,
        history_size::Int = 3,
    ) where {T <: Real}
    n_bands(model) != n_wannier(model) &&
        error("n_bands != n_wann, run instead disentanglement?")

    fg! = get_fg!_maxloc(p, model)

    Ωⁱ = omega(p, model.kstencil, model.overlaps, model.gauges)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    kManif = Optim.Stiefel_SVD()
    Manif = Optim.PowerManifold(kManif, (n_wannier(model), n_wannier(model)), (n_kpoints(model),))

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    Uinit = [
        model.gauges[ik][ib, ic] for ib in 1:(n_bands(model)), ic in 1:(n_wannier(model)),
            ik in 1:length(model.gauges)
    ]

    opt = Optim.optimize(
        Optim.only_fg!(fg!),
        Uinit,
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

    Ωᶠ = omega(p, model.kstencil, model.overlaps, Umin)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    Umin_vec = map(1:(n_kpoints(model))) do ik
        Umin[:, :, ik]
    end

    return Umin_vec
end

localize_isolated_bands(model::Model; kwargs...) =
    localize_isolated_bands(SpreadPenalty(), model; kwargs...)

max_localize(p::AbstractPenalty, model::Model{T}; kwargs...) where {T <: Real} =
    localize_isolated_bands(p, model; kwargs...)

max_localize(model::Model; kwargs...) = localize_isolated_bands(model; kwargs...)
