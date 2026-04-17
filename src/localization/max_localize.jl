export max_localize, localize

"""
    max_localize(model; kwargs...)

Maximally localize the Marzari-Vanderbilt spread on an isolated
manifold. `kwargs` are forwarded to [`OptimLBFGS`](@ref).
"""
max_localize(model::Model; kwargs...) =
    solve!(Problem(Variance(), model), OptimLBFGS(; kwargs...))

"""
    localize(model, r0, λ; kwargs...)

Maximally localize `model` with a per-WF center penalty
`Ωc = λ · Σₙ |r_n − r0[n]|²`.
"""
localize(model::Model, r0::AbstractVector, λ::Real; kwargs...) =
    solve!(Problem(CenteredVariance(collect(Vec3{Float64}, r0), Float64(λ)), model), OptimLBFGS(; kwargs...))
