export localize

"""
    localize(model; kwargs...)
    localize(sm::SpinModel; λs=1.0, kwargs...)
    localize(obj::Objective, model; kwargs...)
    localize(obj::Objective, model, layout::Layout; kwargs...)

Run localization against `model` (or `SpinModel` `sm`). The no-objective
forms default to [`Variance`](@ref) on a `Model` and [`CoOptVariance`](@ref)
on a `SpinModel`; pass a concrete [`Objective`](@ref) to use
[`CenteredVariance`](@ref) / [`CenteredCoOptVariance`](@ref) / etc. The
four-argument form lets the caller pick the packing [`Layout`](@ref)
explicitly (e.g. `WLayout()` for opt-rotate-style single-W optimization).
All `kwargs` are forwarded to [`OptimLBFGS`](@ref).
"""
localize(model::Model; kwargs...) = localize(Variance(), model; kwargs...)
localize(sm::SpinModel; λs::Real = 1.0, kwargs...) =
    localize(CoOptVariance(float(λs)), sm; kwargs...)

localize(obj::Objective, model; kwargs...) =
    solve!(Problem(obj, model), OptimLBFGS(; kwargs...))

localize(obj::Objective, model, layout::Layout; kwargs...) =
    solve!(Problem(obj, model, layout), OptimLBFGS(; kwargs...))
