export constrain_center_coopt

function omega(
        terms::Tuple,
        model::SpinModel,
        Uup::AbstractArray{T, 3},
        Udn::AbstractArray{T, 3},
        λs::R,
    ) where {T <: Complex, R <: Real}
    center_term = _find_center_term(terms)
    isnothing(center_term) && error("CenterConstraintTerm is required for constrained-center magnetic omega")
    up = omega_center(omega(model.up, Uup); r₀ = center_term.r0, λ = center_term.λ)
    dn = omega_center(omega(model.dn, Udn); r₀ = center_term.r0, λ = center_term.λ)
    M = overlap_updn(model, Uup, Udn)
    Ωupdn = omega_updn(M)
    Ωt = up.Ωt + dn.Ωt + λs * Ωupdn
    return SpreadMag(up, dn, Ωupdn, Ωt, M, λs)
end

function omega(terms::Tuple, model::SpinModel{T}, λs::T) where {T <: Real}
    return omega(terms, model, model.up.gauges, model.dn.gauges, λs)
end

"""
    get_fg!_disentangle(p, model::SpinModel, λs)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
function get_fg!_disentangle(
        terms::Tuple, model::SpinModel{T}, λs::T
    ) where {T <: Real}
    center = _find_center_term(terms)
    isnothing(center) &&
        error("CenterConstraintTerm is required for constrained-center SpinModel fg!")
    obj = CenteredCoOptVariance(center.r0, center.λ, float(λs))
    return _make_optim_fg!(Problem(obj, model))
end

"""
    constrain_center_coopt(sm, r0, λ; λs=1.0, kwargs...)

Co-optimize the spin-up and spin-down Wannier gauges of a [`SpinModel`](@ref)
with both the ↑↓ overlap coupling (weighted by `λs`) and a per-WF center
penalty `λ · Σₙ |r_n − r0[n]|²` on each channel.
"""
function constrain_center_coopt(
        sm::SpinModel, r0::AbstractVector, λ::Real; λs::Real = 1.0, kwargs...
    )
    obj = CenteredCoOptVariance(collect(Vec3{Float64}, r0), Float64(λ), Float64(λs))
    return solve!(Problem(obj, sm), OptimLBFGS(; kwargs...))
end

# Back-compat: accept the legacy tuple-terms signature for constrain_center coopt.
function disentangle(
        terms::Tuple, sm::SpinModel, λs::Real = 1.0; kwargs...
    )
    center_term = _find_center_term(terms)
    isnothing(center_term) &&
        error("constrain_center disentangle(terms, sm, λs) requires a CenterConstraintTerm")
    return constrain_center_coopt(sm, center_term.r0, center_term.λ; λs = λs, kwargs...)
end

