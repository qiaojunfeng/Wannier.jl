"""
Model for spin polarized system with constraint.

Traditionally, we run two independent Wannierizations for spin up and spin down.
Here we add a constraint to maximally overlap the spin-up and spin-down WFs,
so that they map one-by-one to each other.

The constructor enforces that `up` and `dn` describe the same real-space /
reciprocal-space system (lattice, atoms, kstencil); the `M` field stores the
Bloch-basis ↑↓ overlap `< u_nk^↑ | u_mk^↓ >` read from the coupling amn file.
"""
struct SpinModel{T <: Real}
    up::Model{T}
    dn::Model{T}

    "< u_nk^↑ | u_mk^↓ >, size: (n_bands, n_bands, n_kpts)"
    M::Array{Complex{T}, 3}

    function SpinModel{T}(up::Model{T}, dn::Model{T}, M::AbstractArray{<:Complex, 3}) where {T <: Real}
        isapprox(up.lattice, dn.lattice) || error("SpinModel: up/dn lattice mismatch")
        (length(up.atom_positions) == length(dn.atom_positions) &&
            all(isapprox(a, b) for (a, b) in zip(up.atom_positions, dn.atom_positions))) ||
            error("SpinModel: up/dn atom_positions mismatch")
        up.atom_labels == dn.atom_labels || error("SpinModel: up/dn atom_labels mismatch")
        isapprox(up.kstencil, dn.kstencil) || error("SpinModel: up/dn kstencil mismatch")
        size(M, 3) == n_kpoints(up) || error("SpinModel: M has wrong number of kpoints")
        size(M, 1) == size(M, 2) == n_bands(up) || error("SpinModel: M has wrong n_bands")
        return new{T}(up, dn, Array{Complex{T}, 3}(M))
    end
end

function SpinModel(up::Model{T}, dn::Model{T}, M::AbstractArray{<:Complex, 3}) where {T <: Real}
    return SpinModel{T}(up, dn, M)
end

function Base.show(io::IO, ::MIME"text/plain", model::SpinModel)
    println(io, "spin up:")
    show(io, model.up)
    println(io, "\n")

    println(io, "spin down:")
    return show(io, model.dn)
end

# Accessors delegate to `up` — constructor invariants guarantee they match `dn`.
n_atoms(m::SpinModel) = n_atoms(m.up)
n_kpoints(m::SpinModel) = n_kpoints(m.up)
n_bvectors(m::SpinModel) = n_bvectors(m.up)
n_bands(m::SpinModel) = n_bands(m.up)
n_wannier(m::SpinModel) = n_wannier(m.up)
CrystalBase.real_lattice(m::SpinModel) = real_lattice(m.up)
CrystalBase.reciprocal_lattice(m::SpinModel) = reciprocal_lattice(m.up)
kpoints(m::SpinModel) = kpoints(m.up)
kgrid_size(m::SpinModel) = kgrid_size(m.up)
kpb_k(m::SpinModel) = kpb_k(m.up)
kpb_G(m::SpinModel) = kpb_G(m.up)
bweights(m::SpinModel) = bweights(m.up)

struct SpreadMag{T <: Real, S <: AbstractSpread} <: AbstractSpread
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

function omega(model::SpinModel, Uup, Udn, λ::Real)
    up = omega(model.up, Uup)
    dn = omega(model.dn, Udn)
    M = overlap_updn(model, Uup, Udn)
    Ωupdn = omega_updn(M)
    Ωt = up.Ω + dn.Ω + λ * Ωupdn
    return SpreadMag(up, dn, Ωupdn, Ωt, M, λ)
end

function omega(model::SpinModel, λ::Real)
    return omega(model, model.up.gauges, model.dn.gauges, λ)
end

function Base.show(io::IO, M::MIME"text/plain", Ω::SpreadMag)
    @info "spin up:"
    show(io, M, Ω.up)
    println(io, "\n")

    @info "spin down:"
    show(io, M, Ω.dn)
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
    return @printf(io, "   Ωt  = %11.5f\n", Ω.Ωt)
end

"""
    overlap_updn(M, Uup, Udn)

Compute the overlap between up and down WFs.

Actually N - Ω↑↓, according to QPPM Eq. 8, where N = n_wann.

# Arguments
- `M`: the `SpinModel.M` matrices, `(n_bands, n_bands, n_kpts)`
- `Uup`: the up gauge matrices, `(n_bands, n_wann, n_kpts)`
- `Udn`: the down gauge matrices, `(n_bands, n_wann, n_kpts)`
"""
function overlap_updn(
        M::AbstractArray{<:Complex, 3},
        Uup::AbstractArray{<:Complex, 3},
        Udn::AbstractArray{<:Complex, 3},
    )
    n_bands, n_wann, n_kpts = size(Uup)
    Mᵂ = zeros(eltype(M), n_wann, n_wann)

    for ik in 1:n_kpts
        Mᵂ .+= view(Uup, :, :, ik)' * view(M, :, :, ik) * view(Udn, :, :, ik)
    end

    return abs2.(Mᵂ) ./ n_kpts^2
end

function overlap_updn(model::SpinModel, Uup, Udn)
    return overlap_updn(model.M, Uup, Udn)
end

function overlap_updn(model::SpinModel)
    return overlap_updn(model, model.up.gauges, model.dn.gauges)
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

function omega_updn(model::SpinModel, Uup, Udn)
    return omega_updn(overlap_updn(model, Uup, Udn))
end

function omega_updn(model::SpinModel)
    return omega_updn(overlap_updn(model))
end

@doc raw"""
    overlap_updn_grad(model::SpinModel, Uup, Udn)

Compute gradients of [`overlap_updn`](@ref overlap_updn).

``\frac{d \Omega}{d U^{\uparrow}}`` and ``\frac{d \Omega}{d U^{\downarrow}}``.

TODO: this is actually the gradient of Tr[overlap_updn]

# Arguments
- `M`: the `SpinModel.M` matrices, `(n_bands, n_bands, n_kpts)`
- `Uup`: the up gauge matrices, `(n_bands, n_wann, n_kpts)`
- `Udn`: the down gauge matrices, `(n_bands, n_wann, n_kpts)`
"""
function overlap_updn_grad(
        M::AbstractArray{<:Complex, 3},
        Uup::AbstractArray{<:Complex, 3},
        Udn::AbstractArray{<:Complex, 3},
    )
    n_bands, n_wann, n_kpts = size(Uup)

    T = eltype(Uup)
    Mᵂ = zeros(T, n_wann, n_wann)
    GUup = zeros(T, n_bands, n_wann, n_kpts)
    GUdn = zeros(T, n_bands, n_wann, n_kpts)

    for ik in 1:n_kpts
        Uupk = view(Uup, :, :, ik)
        Udnk = view(Udn, :, :, ik)
        Mk = view(M, :, :, ik)
        MUdn = Mk * Udnk
        Mᵂ += Uupk' * MUdn

        GUup[:, :, ik] .= MUdn
        GUdn[:, :, ik] .= Mk' * Uupk
    end

    diagM = diagm(diag(Mᵂ))
    for ik in 1:n_kpts
        view(GUup, :, :, ik) .= view(GUup, :, :, ik) * diagM'
        view(GUdn, :, :, ik) .= view(GUdn, :, :, ik) * diagM
    end
    return GUup ./ n_kpts^2, GUdn ./ n_kpts^2
end

function overlap_updn_grad(model::SpinModel, Uup, Udn)
    return overlap_updn_grad(model.M, Uup, Udn)
end

@doc raw"""
    overlap_updn_grad(model::SpinModel, Xup, Yup, Xdn, Ydn)

Compute gradient of [`overlap_updn`](@ref overlap_updn).

``\frac{d \Omega}{d X^{\uparrow}}``, ``\frac{d \Omega}{d Y^{\uparrow}}``,
``\frac{d \Omega}{d X^{\downarrow}}``, ``\frac{d \Omega}{d Y^{\downarrow}}``.
"""
function overlap_updn_grad(model::SpinModel, Xup, Yup, Xdn, Ydn)
    Uup = X_Y_to_U(Xup, Yup)
    Udn = X_Y_to_U(Xdn, Ydn)
    GUup, GUdn = overlap_updn_grad(model, Uup, Udn)

    GXup, GYup = GU_to_GX_GY(GUup, Xup, Yup, model.up.frozen_bands)
    GXdn, GYdn = GU_to_GX_GY(GUdn, Xdn, Ydn, model.dn.frozen_bands)

    return GXup, GYup, GXdn, GYdn
end

function omega_updn_grad(model::SpinModel, Uup, Udn)
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

function omega_updn_grad(model::SpinModel, Xup, Yup, Xdn, Ydn)
    GXup, GYup, GXdn, GYdn = overlap_updn_grad(model, Xup, Yup, Xdn, Ydn)
    # see previous function for the reason of the factor of -2
    return -2 * GXup, -2 * GYup, -2 * GXdn, -2 * GYdn
end

"""
    get_fg!_disentangle(model::SpinModel)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_disentangle(model::SpinModel, λ::Real = 1.0)
    problem = LocalizationProblem((VarianceTerm(),), model, :mag_disentangle; lambda = λ)
    return build_fg!(problem)
end

"""
    disentangle(model; f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Run disentangle on a `SpinModel`.

# Arguments
- `model`: SpinModel
- `λ`: Lagrange multiplier of the ↑↓ overlap term

# Keyword arguments
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function disentangle(
        model::SpinModel{T},
        λ::T = 1.0;
        f_tol::T = 1.0e-7,
        g_tol::T = 1.0e-5,
        max_iter::Int = 200,
        history_size::Int = 3,
    ) where {T <: Real}
    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)

    @assert n_bands(model.dn) == nb
    @assert n_wannier(model.dn) == nw
    @assert n_kpoints(model.dn) == nk

    XYk_up_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw)
    )
    XYk_dn_Manif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (nw, nw), (nb, nw)
    )
    n_inner = nw^2 + nb * nw
    XYkManif = Optim.ProductManifold(XYk_up_Manif, XYk_dn_Manif, (n_inner,), (n_inner,))
    XYManif = Optim.PowerManifold(XYkManif, (2 * n_inner,), (nk,))

    Xup0, Yup0 = U_to_X_Y(model.up.gauges, model.up.frozen_bands)
    Xdn0, Ydn0 = U_to_X_Y(model.dn.gauges, model.dn.frozen_bands)
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

    ls = Optim.HagerZhang()
    meth = Optim.LBFGS

    opt = Optim.optimize(
        f,
        g!,
        XY0,
        meth(; manifold = XYManif, linesearch = ls, m = history_size),
        Optim.Options(;
            show_trace = true,
            iterations = max_iter,
            f_tol = f_tol,
            g_tol = g_tol,
            allow_f_increases = true,
        ),
    )
    display(opt)

    XYmin = Optim.minimizer(opt)

    XYupmin = XYmin[1:n_inner, :]
    XYdnmin = XYmin[(n_inner + 1):end, :]
    Xupmin, Yupmin = XY_to_X_Y(XYupmin, nb, nw)
    Xdnmin, Ydnmin = XY_to_X_Y(XYdnmin, nb, nw)
    Uupmin = X_Y_to_U(Xupmin, Yupmin)
    Udnmin = X_Y_to_U(Xdnmin, Ydnmin)

    @info "Final spread"
    Ω = omega(model, Uupmin, Udnmin, λ)
    show(Ω)

    return Uupmin, Udnmin
end
