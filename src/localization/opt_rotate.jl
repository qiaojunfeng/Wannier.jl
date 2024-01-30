using LinearAlgebra
using Optim: Optim
import NLSolversBase: OnceDifferentiable

export opt_rotate, merge_gauge

"""
    get_fg!_rotate(model::Model)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_rotate(model::Model)
    function f(W)
        U = merge_gauge(model.gauges, W)
        return omega(model.kstencil, model.overlaps, U).Ω
    end

    function g!(G, W)
        n_wann = size(W, 1)
        M = model.overlaps
        n_bvecs = length(M[1])
        n_kpts = n_kpoints(model)

        bvectors = model.kstencil
        kpb_k = bvectors.kpb_k
        kpb_G = bvectors.kpb_G
        wb = bvectors.bweights
        recip_lattice = bvectors.recip_lattice
        kpoints = bvectors.kpoints

        fill!(G, 0.0)
        T = zeros(eltype(W), n_wann, n_wann)
        MWᵏᵇ = zeros(eltype(W), n_wann, n_wann)
        Nᵏᵇ = zeros(eltype(W), n_wann, n_wann)

        UW = merge_gauge(model.gauges, W)
        r = center(bvectors, M, UW)
        # actually I can just call this, equivalent to the for loop below.
        # G_U = omega_grad(bvectors, M, UW, r)
        # # sum w.r.t. kpoints
        # G .= dropdims(sum(G_U, dims = 3); dims = 3)
        # return nothing

        for ik in 1:n_kpts
            for ib in 1:n_bvecs
                ikpb = kpb_k[ik][ib]

                # need to use UW[:, :, ik] instead of W, if model.U is not identity
                MWᵏᵇ .= M[ik][ib] * W
                Nᵏᵇ .= W' * MWᵏᵇ
                b = recip_lattice * (kpoints[ikpb] + kpb_G[ik][ib] - kpoints[ik])

                q = imaglog.(diag(Nᵏᵇ))
                for iw in 1:n_wann
                    q[iw] += r[iw] ⋅ b
                end

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
    opt_rotate(model; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Maximally localize spread functional w.r.t. single unitary matrix `W`.

# Arguments
- `model`: model

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function opt_rotate(
    model::Model{T}; f_tol::T=1e-7, g_tol::T=1e-5, max_iter::Int=200, history_size::Int=3
) where {T<:Real}
    n_wann = n_wannier(model)
    n_bands(model) == n_wann || error("n_bands != n_wann, run instead disentanglement?")

    Ωⁱ = omega(model)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    # make sure U is identity matrices
    # however I shouldn't modify the original model, deepcopy it
    # note I cannot use `rotate_gauge(model, model.U)` because the rotated Hamiltonian
    # might not be diagonal, however in this case I don't care about eigenvalues
    model2 = deepcopy(model)
    model2.overlaps .= transform_gauge(
        model2.overlaps, model2.kstencil.kpb_k, model2.gauges
    )
    model2.gauges .= identity_gauge(eltype(model2.gauges[1]), n_kpoints(model2), n_wann)

    wManif = Optim.Stiefel_SVD()

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    W0 = Matrix{eltype(model2.gauges[1])}(I, n_wann, n_wann)

    f, g! = get_fg!_rotate(model2)

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

    # model2.U is actually identity
    U = merge_gauge(model2.gauges, Wmin)
    Ωᶠ = omega(model2, U)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Wmin
end

"""
    merge_gauge(U::Array{T,3}, W::Matrix{T}) where {T<:Complex}

Rotate the `U` matrices at each kpoint by the same `W` matrix.

``\\forall \\bm{k}``, ``U_{\\bm{k}} W``

Useful once we have the optimal rotation matrix `W`, then update the initial
`U` matrices by rotating them by `W`.
"""
function merge_gauge(U::Vector, W::Matrix{T}) where {T<:Complex}
    n_bands, n_wann = size(U[1])

    size(W) != (n_wann, n_wann) && error("W must be a n_wann x n_wann matrix")
    return map(u -> u * W, U)
end

function merge_gauge(U::Array{T,3}, W::Matrix{T}) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(U)
    size(W) != (n_wann, n_wann) && error("W must be a n_wann x n_wann matrix")

    U1 = similar(U)

    for ik in 1:n_kpts
        U1[:, :, ik] .= U[:, :, ik] * W
    end

    return U1
end
