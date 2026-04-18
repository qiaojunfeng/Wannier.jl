export ParallelTransport

"""
    ParallelTransport(; use_U=false, log_interp=false)

Closed-form gauge construction by parallel transport along a MP-like
k-grid plus obstruction pullbacks (GLS2019). Produces a smooth gauge
from an initial one *without* minimizing a scalar functional; there is
no `Problem`, no `fg!`, no solver. `localize(ParallelTransport(), model)`
routes to the underlying [`parallel_transport`](@ref) algorithm.

# Keyword arguments
- `use_U`: seed from the current `model.gauges` instead of the identity.
- `log_interp`: use logarithmic interpolation for the line/surface
  obstruction pullbacks.
"""
Base.@kwdef struct ParallelTransport
    use_U::Bool = false
    log_interp::Bool = false
end
