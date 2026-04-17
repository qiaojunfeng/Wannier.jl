export constrain_center_coopt

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
