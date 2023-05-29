using Optim: Optim

export max_localize 

"""
    get_fg!_maxloc(model::Model)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_maxloc(p::AbstractPenalty, model::Model)
    cache = Cache(model)

    function fg!(F, G, U)
        compute_MUᵏᵇ_Nᵏᵇ!(cache, model.bvectors, model.M, U)
        
        if G !== nothing
            cache.G = G
            omega_grad!(p, cache, model.bvectors, model.M)
        end
        if F !== nothing
            return omega!(p, cache, model.bvectors, model.M).Ω
        end
    end

    return fg!
end

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
function max_localize(p::AbstractPenalty,
    model::Model{T}; f_tol::T=1e-7, g_tol::T=1e-5, max_iter::Int=200, history_size::Int=3
) where {T<:Real}
    model.n_bands != model.n_wann &&
        error("n_bands != n_wann, run instead disentanglement?")

    fg! = get_fg!_maxloc(p, model)

    Ωⁱ = omega(p, model.bvectors, model.M, model.U)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    kManif = Optim.Stiefel_SVD()
    Manif = Optim.PowerManifold(kManif, (model.n_wann, model.n_wann), (model.n_kpts,))

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    Uinit = [model.U[ik][ib, ic] for ib=1:size(model.U[1],1), ic = 1:size(model.U[1],2), ik = 1:length(model.U)]

    opt = Optim.optimize(
        Optim.only_fg!(fg!),
        Uinit,
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

    Umin = Optim.minimizer(opt)

    Ωᶠ = omega(p, model.bvectors, model.M, Umin)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Umin
end
