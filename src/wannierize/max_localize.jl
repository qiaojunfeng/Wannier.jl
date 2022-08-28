using Optim: Optim

export max_localize

"""
    get_fg!_maxloc(model::Model)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_maxloc(model::Model)
    f(A) = omega(model.bvectors, model.M, A).Ω

    g!(G, A) = begin
        G .= omega_grad(model.bvectors, model.M, A)
        nothing
    end

    return f, g!
end

"""
    max_localize(model; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=20)

Maximally localize spread functional w.r.t. all kpoints on a unitary matrix manifold.

# Arguments
- `model`: model

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function max_localize(
    model::Model{T}; f_tol::T=1e-7, g_tol::T=1e-5, max_iter::Int=200, history_size::Int=20
) where {T<:Real}
    model.n_bands != model.n_wann &&
        error("n_bands != n_wann, run instead disentanglement?")

    f, g! = get_fg!_maxloc(model)

    Ωⁱ = omega(model.bvectors, model.M, model.A)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    kManif = Optim.Stiefel_SVD()
    Manif = Optim.PowerManifold(kManif, (model.n_wann, model.n_wann), (model.n_kpts,))

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    Ainit = deepcopy(model.A)

    opt = Optim.optimize(
        f,
        g!,
        Ainit,
        meth(; manifold=Manif, linesearch=ls, m=history_size),
        Optim.Options(;
            show_trace=true,
            iterations=max_iter,
            f_tol=f_tol,
            g_tol=g_tol,
            allow_f_increases=true,
        ),
    )
    display(opt)

    Amin = Optim.minimizer(opt)

    Ωᶠ = omega(model.bvectors, model.M, Amin)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Amin
end
