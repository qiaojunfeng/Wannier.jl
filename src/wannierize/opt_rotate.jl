using LinearAlgebra
using Optim: Optim
import NLSolversBase: OnceDifferentiable

function get_fg!_rotate(model::Model)
    function f(W)
        A = rotate_amn(model.A, W)
        return omega(model.bvectors, model.M, A).Ω
    end

    function g!(G, W)
        n_wann = size(W, 1)
        M = model.M
        n_bvecs = size(M, 3)
        n_kpts = model.n_kpts

        bvectors = model.bvectors
        kpb_k = bvectors.kpb_k
        kpb_b = bvectors.kpb_b
        wb = bvectors.weights
        recip_lattice = bvectors.recip_lattice
        kpoints = bvectors.kpoints

        fill!(G, 0.0)
        T = zeros(eltype(W), n_wann, n_wann)
        b = zeros(eltype(real(W)), 3)
        MWᵏᵇ = zeros(eltype(W), n_wann, n_wann)
        Nᵏᵇ = zeros(eltype(W), n_wann, n_wann)

        AW = rotate_amn(model.A, W)
        r = center(bvectors, M, AW)
        # actually I can just call this, equivalent to the for loop below.
        # G_A = omega_grad(bvectors, M, AW, r)
        # # sum w.r.t. kpoints
        # G .= dropdims(sum(G_A, dims = 3); dims = 3)
        # return nothing

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                ikpb = kpb_k[ib, ik]

                MWᵏᵇ .= overlap(M, kpb_k, ik, ikpb) * W
                Nᵏᵇ .= W' * MWᵏᵇ
                b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])

                q = imaglog.(diag(Nᵏᵇ)) + r' * b

                for n in 1:n_wann
                    # error if division by zero. Should not happen if the initial gauge is not too bad
                    if abs(Nᵏᵇ[n, n]) < 1e-10
                        error("Nᵏᵇ too small! $ik -> $ikpb, $Nᵏᵇ")
                    end

                    t = -conj(Nᵏᵇ[n, n]) - im * q[n] / Nᵏᵇ[n, n]

                    for m in 1:n_wann
                        T[m, n] = t * MWᵏᵇ[m, n]
                    end
                end

                G .+= 4 * wb[ib] * T
            end
        end

        G ./= n_kpts

        return nothing
    end

    return f, g!
end

"""
Maximally localize spread functional w.r.t. single unitary matrix W.
"""
function opt_rotate(
    model::Model{T}; f_tol::T=1e-7, g_tol::T=1e-5, max_iter::Int=200, history_size::Int=20
) where {T<:Real}
    n_wann = model.n_wann

    model.n_bands != n_wann && error("n_bands != n_wann, run instead disentanglement?")

    f, g! = get_fg!_rotate(model)

    Ωⁱ = omega(model.bvectors, model.M, model.A)
    @info "Initial spread"
    pprint(Ωⁱ)

    # make sure A is identity matrices
    model.M .= rotate_mmn(model.M, model.bvectors.kpb_k, model.A)
    model.A .= ones_amn(eltype(model.A), model.n_wann, model.n_kpts)

    wManif = Optim.Stiefel_SVD()

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    W0 = Matrix{eltype(model.A)}(I, n_wann, n_wann)

    # Comment g! to use finite differences to compute the gradient
    opt = Optim.optimize(
        f,
        g!,
        W0,
        meth(; manifold=wManif, linesearch=ls, m=history_size),
        # autodiff=:forward,
        Optim.Options(;
            show_trace=true,
            iterations=max_iter,
            f_tol=f_tol,
            g_tol=g_tol,
            allow_f_increases=true,
        ),
    )
    display(opt)

    Wmin = Optim.minimizer(opt)

    A = rotate_amn(model.A, Wmin)
    Ωᶠ = omega(model.bvectors, model.M, A)
    @info "Final spread"
    pprint(Ωᶠ)

    return Wmin
end

function rotate_amn(A::Array{T,3}, W::Matrix{T}) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(A)
    size(W) != (n_wann, n_wann) && error("W must be a n_wann x n_wann matrix")

    A1 = similar(A)

    for ik in 1:n_kpts
        A1[:, :, ik] .= A[:, :, ik] * W
    end

    return A1
end
