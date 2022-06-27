import LinearAlgebra as LA
import Optim
import NLSolversBase: OnceDifferentiable


function get_fg!_rotate(model::Model)

    function f(W)
        A = similar(model.A)

        for i = 1:model.n_kpts
            A[:, :, i] .= model.A[:, :, i] * W
        end

        omega(model.bvectors, model.M, A).Ω
    end

    function g!(G, W)
        A = similar(model.A)

        for i = 1:model.n_kpts
            A[:, :, i] .= model.A[:, :, i] * W
        end

        GA = omega_grad(model.bvectors, model.M, A)

        # sum w.r.t. kpoints
        G .= dropdims(sum(GA, dims = 3); dims = 3)

        nothing
    end

    f, g!
end


"""
Maximally localize spread functional w.r.t. single unitary matrix W.
"""
function opt_rotate(
    model::Model{T};
    f_tol::T = 1e-10,
    g_tol::T = 1e-8,
    max_iter::Int = 1000,
    history_size::Int = 20,
) where {T<:Real}
    n_wann = model.n_wann

    model.n_bands != n_wann && error("n_bands != n_wann, run instead disentanglement?")

    f, g! = get_fg!_rotate(model)

    Ωⁱ = f(model.A)
    @info "Initial total spread" Ω = round(Ωⁱ.Ω; digits = 5)
    @info "Initial spread:" ω = round.(Ωⁱ.ω'; digits = 5)
    @info "Initial centers:" r = round.(Ωⁱ.r; digits = 5)

    wManif = Optim.Stiefel_SVD()

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    W0 = Matrix{eltype(model.A)}(LA.I, n_wann, n_wann)

    # use Finite differences to compute the gradient
    opt = Optim.optimize(
        f,
        #g!,
        W0,
        meth(manifold = wManif, linesearch = ls, m = history_size),
        # autodiff=:forward,
        Optim.Options(
            show_trace = true,
            iterations = max_iter,
            f_tol = f_tol,
            g_tol = g_tol,
            allow_f_increases = true,
        ),
    )
    display(opt)

    Wmin = Optim.minimizer(opt)

    # Wmin = [  0.396394+0.304749im  -0.0652139-0.49573im     -0.382338+0.322204im  -0.420344-0.27076im
    # -0.254709-0.410673im    0.702612-0.0688952im  -0.0219307+0.433408im   -0.26895+0.086048im
    #  0.526111+0.233063im    0.429235+0.0678782im    0.624165-0.22825im   -0.187253-0.0574045im
    #  0.428172-0.04504im     0.182706-0.17134im    -0.0799552+0.337054im     0.7812+0.147068im]

    Amin = similar(model.A)
    for i = 1:model.n_kpts
        Amin[:, :, i] .= model.A[:, :, i] * Wmin
    end

    Ωᶠ = f(Amin)
    @info "Final total spread" Ω = round(Ωᶠ.Ω; digits = 5)
    @info "Final spread:" ω = round.(Ωᶠ.ω'; digits = 5)
    @info "Final centers:" r = round.(Ωᶠ.r; digits = 5)

    Amin
end
