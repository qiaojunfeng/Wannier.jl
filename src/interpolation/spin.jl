export SpinRspace, SpinKspace

"""
    $(TYPEDEF)

A struct representing tight-binding spin operator in R-space.

# Fields
$(FIELDS)

!!! note

    This is defined in the same way as [`PositionRspace`](@ref), however,
    since the spin operator and position operator transform differently
    under gauge transformation, we define them separately so that they can
    be dispatched differently.
"""
struct SpinRspace{M<:AbstractMatrix} <: AbstractOperatorRspace
    """the R-space domain (or called R-vectors) on which the operator is defined"""
    domain::BareRspaceDomain

    """The tight-binding operator defined on the domain."""
    operator::Vector{M}
end

"""
    $(TYPEDEF)

A struct representing tight-binding spin operator on a kpoint grid.

# Fields
$(FIELDS)
"""
struct SpinKspace{K<:AbstractKpointContainer,M<:AbstractMatrix} <: AbstractOperatorKspace{K}
    """a [`KpointGrid`](@ref) or [`KpointList`](@ref) on which the operator is defined"""
    domain::K

    """the tight-binding operator defined on the domain"""
    operator::Vector{M}
end
