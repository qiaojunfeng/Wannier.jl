"""
Model for spin polarized system with constraint.

Traditionally, we run two independent Wannierizations for spin up and spin down.
Here we add a constraint to maximally overlap the spin-up and spin-down WFs,
so that they map one-by-one to each other.
"""
struct MagModel{T<:Real}
    # submodel for spin up
    up::Model{T}

    # submodel for spin down
    dn::Model{T}

    # < u_nk^↑ | u_mk^↓ >, size: (n_bands, n_bands, n_kpts)
    M::Vector{Matrix{Complex{T}}}
end

function Base.show(io::IO, model::MagModel)
    # TODO not sure why but this breaks REPL
    # every time I run any command, it will show this info message?
    # @info "spin up:"
    println(io, "spin up:")
    show(io, model.up)
    println(io, "\n")

    # @info "spin down:"
    println(io, "spin down:")
    return show(io, model.dn)
end

# function MagModel(up::Model{T}, dn::Model{T}) where {T<:Real}
#     @assert up.n_bands == dn.n_bands
#     @assert up.n_wann == dn.n_wann
#     @assert up.n_kpts == dn.n_kpts
#     @assert up.n_bvecs == dn.n_bvecs
#     @assert up.kgrid == dn.kgrid
#     @assert up.kpoints == dn.kpoints
#     @assert up.bvectors == dn.bvectors
#     # @assert up.frozen_bands == dn.frozen_bands

#     MagModel(up, dn)
# end

struct SpreadMag{T<:Real,S<:AbstractSpread} <: AbstractSpread
    # up spread
    up::S

    # dn spread
    dn::S

    # unit Å²
    Ωupdn::T

    # Total spread Ωt = up.Ω + dn.Ω + λ * Ω↑↓
    Ωt::T

    # overlap matrix between up and down WFs, unit Å², size = (n_wann, n_wann)
    M::Matrix{T}

    # λ
    λ::T
end

function omega(model::MagModel, Uup, Udn, λ::Real) 
    up = omega(model.up, Uup)
    dn = omega(model.dn, Udn)
    M = overlap_updn(model, Uup, Udn)
    Ωupdn = omega_updn(M)
    Ωt = up.Ω + dn.Ω + λ * Ωupdn
    return SpreadMag(up, dn, Ωupdn, Ωt, M, λ)
end

function omega(model::MagModel, λ::Real)
    return omega(model, model.up.U, model.dn.U, λ)
end

function Base.show(io::IO, Ω::SpreadMag)
    @info "spin up:"
    show(io, Ω.up)
    println(io, "\n")

    @info "spin down:"
    show(io, Ω.dn)
    println(io, "\n")

    n_wann = size(Ω.M, 1)
    @info "overlap between up and down WFs:"
    @printf(io, "  WF     <↑|↓>/Å²\n")
    for i in 1:n_wann
        @printf(io, "%4d %11.5f\n", i, Ω.M[i, i])
    end
    println(io, "")

    @info "Sum spread: Ωt = Ω↑ + Ω↓ + λ * Ω↑↓"
    @printf(io, "   Ω↑  = %11.5f\n", Ω.up.Ω)
    @printf(io, "   Ω↓  = %11.5f\n", Ω.dn.Ω)
    @printf(io, "   Ω↑↓ = %11.5f\n", Ω.Ωupdn)
    @printf(io, "   Ωt  = %11.5f\n", Ω.Ωt)
end

"""
    overlap_updn(up::Model{T}, dn::Model{T}, M::Array{Complex{T},3}) where {T<:Real}

Compute the overlap between up and down WFs.

Actually N - Ω↑↓, according to QPPM Eq. 8, where N = n_wann.

# Arguments
- `M`: the `MagModel.M` matrices, size (n_bands, n_bands, n_kpts)
- `Uup`: the up gauge matrices, size: (n_bands, n_wann, n_kpts)
- `Udn`: the down gauge matrices, size: (n_bands, n_wann, n_kpts)
"""
function overlap_updn(M::Vector, Uup::Vector, Udn::Vector)
    n_bands, n_wann = size(Uup[1])
    n_kpts = length(Uup)
    Mᵂ = zeros(eltype(M[1]), n_wann, n_wann)

    for ik in 1:n_kpts
        Uupk = Uup[ik]
        Udnk = Udn[ik]
        Mk = M[ik]
        Mᵂ += Uupk' * Mk * Udnk
    end

    return abs2.(Mᵂ) / n_kpts^2
end

function overlap_updn(model::MagModel, Uup, Udn)
    return overlap_updn(model.M, Uup, Udn)
end

function overlap_updn(model::MagModel)
    return overlap_updn(model, model.up.U, model.dn.U)
end

"""
    omega_updn(M)

Compute QPPM Eq. 8.

# Arguments
- `M`: the overlap matrix between up and down WFs, size: (n_wann, n_wann),
    should be the matrix returned from [`overlap_updn`](@ref).
"""
function omega_updn(M::AbstractMatrix)
    n_wann = size(M, 1)
    # I am using minus sign here because the optimizer is minimizing total spread,
    # thus maximizing the ↑↓ overlap.
    return n_wann - sum(diag(M))
end

function omega_updn(model::MagModel, Uup, Udn)
    return omega_updn(overlap_updn(model, Uup, Udn))
end

function omega_updn(model::MagModel)
    return omega_updn(overlap_updn(model))
end

@doc raw"""
    overlap_updn_grad(model::MagModel, Uup, Udn)

Compute gradients of [`overlap_updn`](@ref overlap_updn).

``\frac{d \Omega}{d U^{\uparrow}}`` and ``\frac{d \Omega}{d U^{\downarrow}}``.

TODO: this is actually the gradient of Tr[overlap_updn]

# Arguments
- `M`: the `MagModel.M` matrices, size (n_bands, n_bands, n_kpts)
- `Uup`: the up gauge matrices, size: (n_bands, n_wann, n_kpts)
- `Udn`: the down gauge matrices, size: (n_bands, n_wann, n_kpts)
"""
function overlap_updn_grad(M::Vector, Uup::Vector, Udn::Vector)
    n_bands, n_wann = size(Uup[1])
    n_kpts = length(Uup)

    T = eltype(Uup[1])
    Mᵂ = zeros(T, n_wann, n_wann)
    GUup = [zeros(T, size(Uup[1])) for i = 1:n_kpts]
    GUdn = [zeros(T, size(Udn[1])) for i = 1:n_kpts]

    for ik in 1:n_kpts
        Uupk = Uup[ik]
        Udnk = Udn[ik]
        Mk = M[ik]
        MUdn = Mk * Udnk
        Mᵂ += Uupk' * MUdn

        GUup[ik] = MUdn
        GUdn[ik] = Mk' * Uupk
    end

    diagM = diagm(diag(Mᵂ))
    for ik in 1:n_kpts
        GUup[ik] .= GUup[ik] * diagM'
        GUdn[ik] .= GUdn[ik] * diagM
    end
    return GUup ./ n_kpts^2, GUdn ./ n_kpts^2
end

function overlap_updn_grad(model::MagModel, Uup, Udn)
    return overlap_updn_grad(model.M, Uup, Udn)
end

@doc raw"""
    overlap_updn_grad(model::MagModel, Xup, Yup, Xdn, Ydn)

Compute gradient of [`overlap_updn`](@ref overlap_updn).

``\frac{d \Omega}{d X^{\uparrow}}``, ``\frac{d \Omega}{d Y^{\uparrow}}``,
``\frac{d \Omega}{d X^{\downarrow}}``, ``\frac{d \Omega}{d Y^{\downarrow}}``.
"""
function overlap_updn_grad(
    model::MagModel,
    Xup,
    Yup,
    Xdn,
    Ydn,
)
    Uup = X_Y_to_U(Xup, Yup)
    Udn = X_Y_to_U(Xdn, Ydn)
    GUup, GUdn = overlap_updn_grad(model, Uup, Udn)

    GXup, GYup = GU_to_GX_GY(GUup, Xup, Yup, model.up.frozen_bands)
    GXdn, GYdn = GU_to_GX_GY(GUdn, Xdn, Ydn, model.dn.frozen_bands)

    return GXup, GYup, GXdn, GYdn
end

function omega_updn_grad(model::MagModel, Uup, Udn)
    GUup, GUdn = overlap_updn_grad(model, Uup, Udn)
    # Note since both Optim.jl (used in the actual minimization) and
    # NLSolversBase.jl (used in test for finite difference check) are adopting
    # the convention of
    #   df(x) = Re<∇f, dx>
    # for the complex differentials, I need to multiply by 2 here, since
    # I am using the convention of
    #   df(x) = 2 Re<∇f, dx>
    # when deriving the gradient by hand, this allows me to reuse existing
    # derivative rules but I need to multiply by 2 here.
    # The minus sign is due to the definition of `omega_updn`.
    return -2 .* GUup, -2 .* GUdn
end

function omega_updn_grad(
    model::MagModel,
    Xup,
    Yup,
    Xdn,
    Ydn,
)
    GXup, GYup, GXdn, GYdn = overlap_updn_grad(model, Xup, Yup, Xdn, Ydn)
    # see previous function for the reason of the factor of -2
    return -2 * GXup, -2 * GYup, -2 * GXdn, -2 * GYdn
end

"""
    get_fg!_disentangle(model::MagModel)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_disentangle(model::MagModel, λ::Real=1.0)
    n_bands = model.up.n_bands
    n_wann = model.up.n_wann
    n_kpts = model.up.n_kpts
    n_inner = n_bands * n_wann + n_wann^2  # size of XY at each k-point

    function f(XY)
        XY = reshape(XY, (2 * n_inner, n_kpts))  # *2 for spin up and down
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, n_bands, n_wann)
        Xdn, Ydn = XY_to_X_Y(XYdn, n_bands, n_wann)
        Ωup = omega(model.up.bvectors, model.up.M, Xup, Yup).Ω
        Ωdn = omega(model.dn.bvectors, model.dn.M, Xdn, Ydn).Ω
        if λ == 0
            Ωupdn = 0
        else
            Ωupdn = omega_updn(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
        end
        return Ωup + Ωdn + λ * Ωupdn
    end

    """size(G) == size(XY)"""
    function g!(G, XY)
        XY = reshape(XY, (2 * n_inner, n_kpts))  # *2 for spin up and down
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, n_bands, n_wann)
        Xdn, Ydn = XY_to_X_Y(XYdn, n_bands, n_wann)
        GXup, GYup = omega_grad(
            model.up.bvectors, model.up.M, Xup, Yup, model.up.frozen_bands
        )
        GXdn, GYdn = omega_grad(
            model.dn.bvectors, model.dn.M, Xdn, Ydn, model.dn.frozen_bands
        )

        # gradient of ↑↓ overlap term
        if λ != 0
            GOXup, GOYup, GOXdn, GOYdn = omega_updn_grad(model, Xup, Yup, Xdn, Ydn)
            GXup += λ * GOXup
            GYup += λ * GOYup
            GXdn += λ * GOXdn
            GYdn += λ * GOYdn
        end

        n = n_wann^2

        for ik in 1:n_kpts
            for (i, v) in enumerate(GXup[ik])
                G[i, ik] = v
            end

            for (i, v) in enumerate(GYup[ik])
                G[n + i, ik] = v
            end

            for (i, v) in enumerate(GXdn[ik])
                G[n_inner + i, ik] = v
            end
            
            for (i, v) in enumerate(GYdn[ik])
                G[n_inner + n + i, ik] = v
            end
        end

        return nothing
    end

    return f, g!
end

"""
    disentangle(model; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Run disentangle on a `MagModel`.

# Arguments
- `model`: MagModel
- `λ`: Lagrange multiplier of the ↑↓ overlap term

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function disentangle(
    model::MagModel{T},
    λ::T=1.0;
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
) where {T<:Real}
    n_bands = model.up.n_bands
    n_wann = model.up.n_wann
    n_kpts = model.up.n_kpts

    @assert model.dn.n_bands == n_bands
    @assert model.dn.n_wann == n_wann
    @assert model.dn.n_kpts == n_kpts

    XYk_up_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (n_wann, n_wann), (n_bands, n_wann)
    )
    XYk_dn_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (n_wann, n_wann), (n_bands, n_wann)
    )
    n_inner = n_wann^2 + n_bands * n_wann
    XYkManif = Optim.ProductManifold(XYk_up_Manif, XYk_dn_Manif, (n_inner,), (n_inner,))
    XYManif = Optim.PowerManifold(XYkManif, (2 * n_inner,), (n_kpts,))

    Xup0, Yup0 = U_to_X_Y(model.up.U, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.U, model.dn.frozen_bands)
    # compact storage
    XYup0 = X_Y_to_XY(Xup0, Yup0)
    XYdn0 = X_Y_to_XY(Xdn0, Ydn0)
    XY0 = vcat(XYup0, XYdn0)

    # We have three storage formats:
    # (X, Y): n_wann * n_wann * n_kpts, n_bands * n_wann * n_kpts
    # U: n_bands * n_wann * n_kpts
    # XY: (n_wann * n_wann + n_bands * n_wann) * n_kpts
    f, g! = get_fg!_disentangle(model, λ)

    @info "Initial spread"
    Ω = omega(model, λ)
    show(Ω)
    println("\n")

    @info "Initial spread (with states freezed)"
    Ω = omega(model, X_Y_to_U(Xup0, Yup0), X_Y_to_U(Xdn0, Ydn0), λ)
    show(Ω)
    println("\n")

    # stepsize_mult = 1
    # step = 0.5/(4*8*p.wb)*(p.N1*p.N2*p.N3)*stepsize_mult
    # ls = LineSearches.Static(step)
    ls = Optim.HagerZhang()
    # ls = LineSearches.BackTracking()

    # meth = Optim.GradientDescent
    # meth = Optim.ConjugateGradient
    meth = Optim.LBFGS

    opt = Optim.optimize(
        f,
        g!,
        XY0,
        meth(; manifold=XYManif, linesearch=ls, m=history_size),
        Optim.Options(;
            show_trace=true,
            iterations=max_iter,
            f_tol=f_tol,
            g_tol=g_tol,
            allow_f_increases=true,
        ),
    )
    display(opt)

    XYmin = Optim.minimizer(opt)

    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, n_bands, n_wann)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, n_bands, n_wann)
    Uupmin = X_Y_to_U(Xupmin, Yupmin)
    Udnmin = X_Y_to_U(Xdnmin, Ydnmin)

    @info "Final spread"
    Ω = omega(model, Uupmin, Udnmin, λ)
    show(Ω)

    return Uupmin, Udnmin
end
