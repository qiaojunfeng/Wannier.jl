export localize

"""
    localize(model; kwargs...)
    localize(sm::SpinModel; λ_spin=1.0, kwargs...)
    localize(sm::SymmetricModel; kwargs...)
    localize(obj::Objective, model; kwargs...)
    localize(obj::Objective, model, layout::Layout; kwargs...)
    localize(pt::ParallelTransport, model)

Run localization against `model` (or `SpinModel` `sm`). Two independent
dispatch paths live here:

- `Objective` methods ([`Variance`](@ref), [`CenteredVariance`](@ref),
  [`CoOptVariance`](@ref), [`CenteredCoOptVariance`](@ref)) minimize a
  scalar spread functional; the call routes through [`Problem`](@ref) +
  [`solve!`](@ref) with [`OptimLBFGS`](@ref). `kwargs` forward to the
  solver. The four-argument form picks the packing [`Layout`](@ref)
  explicitly (e.g. `WLayout()` for opt-rotate-style single-W
  optimization).
- [`ParallelTransport`](@ref) is a closed-form geometric construction —
  no scalar functional, no solver. The call routes through
  [`parallel_transport`](@ref) directly.

The no-method forms default to [`Variance`](@ref) on a `Model` and
[`CoOptVariance`](@ref) on a `SpinModel`. A [`SymmetricModel`](@ref) routes
to `Variance` with the symmetry-constrained [`SymmetricXYLayout`](@ref) and returns
the `(U_fbz, U_ibz)` gauge pair.
"""
localize(model::Model; kwargs...) = localize(Variance(), model; kwargs...)
localize(sm::SpinModel; λ_spin::Real = 1.0, kwargs...) =
    localize(CoOptVariance(float(λ_spin)), sm; kwargs...)

localize(obj::Objective, model; kwargs...) =
    solve!(Problem(obj, model), OptimLBFGS(; kwargs...))

localize(obj::Objective, model, layout::Layout; kwargs...) =
    solve!(Problem(obj, model, layout), OptimLBFGS(; kwargs...))

function localize(pt::ParallelTransport, model::Model; kwargs...)
    U, _ = parallel_transport(model; use_U = pt.use_U, log_interp = pt.log_interp)
    return U
end
