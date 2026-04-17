abstract type AbstractLocalizationTerm end

export VarianceTerm, CenterConstraintTerm

"""
Legacy term marker retained for back-compat with callers (including tests)
that still build `(VarianceTerm(), ...)` tuples. The rewrite routes these
through `Objective` subtypes internally.
"""
struct VarianceTerm <: AbstractLocalizationTerm end

"""Legacy center-penalty term marker (see [`VarianceTerm`](@ref))."""
struct CenterConstraintTerm{T <: Real} <: AbstractLocalizationTerm
    r0::Vector{Vec3{T}}
    λ::T
end

@inline _is_variance_term(::VarianceTerm) = true
@inline _is_variance_term(::AbstractLocalizationTerm) = false
@inline _is_center_term(::CenterConstraintTerm) = true
@inline _is_center_term(::AbstractLocalizationTerm) = false

_has_variance_term(terms) = any(_is_variance_term, terms)
_has_center_term(terms) = any(_is_center_term, terms)

function _find_center_term(terms)
    for term in terms
        if term isa CenterConstraintTerm
            return term
        end
    end
    return nothing
end

"""
    _terms_to_objective(terms)

Translate a legacy `(VarianceTerm(), ...)` tuple into the corresponding
concrete [`Objective`](@ref). `(VarianceTerm(),)` → [`Variance`](@ref);
`(VarianceTerm(), CenterConstraintTerm(r0, λ))` → [`CenteredVariance`](@ref).
"""
function _terms_to_objective(terms)
    _has_variance_term(terms) || error("terms must contain a VarianceTerm")
    center = _find_center_term(terms)
    return isnothing(center) ? Variance() : CenteredVariance(center.r0, center.λ)
end
