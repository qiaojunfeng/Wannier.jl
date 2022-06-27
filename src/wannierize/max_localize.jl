import Optim


function get_fg!_maxloc(model::Model)
    f(A) = omega(model.bvectors, model.M, A).Ω

    g!(G, A) = begin
        G .= omega_grad(model.bvectors, model.M, A)
        nothing
    end

    f, g!
end


"""
Maximally localize spread functional w.r.t. to all kpoints.

On a unitary matrix manifold.
"""
function max_localize(
    model::Model{T};
    f_tol::T = 1e-10,
    g_tol::T = 1e-8,
    max_iter::Int = 1000,
    history_size::Int = 20,
) where {T<:Real}
    model.n_bands != model.n_wann &&
        error("n_bands != n_wann, run instead disentanglement?")

    f, g! = get_fg!_maxloc(model)

    Ωⁱ = f(model.A)
    @info "Initial total spread" Ω = round(Ωⁱ.Ω; digits = 5)
    @info "Initial spread:" ω = round.(Ωⁱ.ω'; digits = 5)
    @info "Initial centers:" r = round.(Ωⁱ.r; digits = 5)

    kManif = Optim.Stiefel_SVD()
    Manif = Optim.PowerManifold(kManif, (model.n_wann, model.n_wann), (model.n_kpts,))

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    Ainit = deepcopy(model.A)

    opt = Optim.optimize(
        f,
        g!,
        Ainit,
        meth(manifold = Manif, linesearch = ls, m = history_size),
        Optim.Options(
            show_trace = true,
            iterations = max_iter,
            f_tol = f_tol,
            g_tol = g_tol,
            allow_f_increases = true,
        ),
    )
    display(opt)

    Amin = Optim.minimizer(opt)

    Ωᶠ = f(Amin)
    @info "Final total spread" Ω = round(Ωᶠ.Ω; digits = 5)
    @info "Final spread:" ω = round.(Ωᶠ.ω'; digits = 5)
    @info "Final centers:" r = round.(Ωᶠ.r; digits = 5)

    Amin
end
