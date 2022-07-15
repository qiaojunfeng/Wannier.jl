using Optim: Optim

function get_fg!_maxloc(model::Model)
    f(A) = omega(model.bvectors, model.M, A).Ω

    g!(G, A) = begin
        G .= omega_grad(model.bvectors, model.M, A)
        nothing
    end

    return f, g!
end

"""
Maximally localize spread functional w.r.t. all kpoints.

On a unitary matrix manifold.
"""
function max_localize(
    model::Model{T}; f_tol::T=1e-7, g_tol::T=1e-5, max_iter::Int=200, history_size::Int=20
) where {T<:Real}
    model.n_bands != model.n_wann &&
        error("n_bands != n_wann, run instead disentanglement?")

    f, g! = get_fg!_maxloc(model)

    Ωⁱ = omega(model.bvectors, model.M, model.A)
    @info "Initial spread"
    pprint(Ωⁱ)

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
    pprint(Ωᶠ)

    return Amin
end
