export TBSpin, SpinInterpolator

"""Construct a tight-binding spin operator in R-space.

!!! note

    This is defined in the same way as [`TBPosition`](@ref). However, note that
    the spin operator and position operator transform differently under gauge
    transformation, see [`SpinInterpolator`](@ref) for details.
"""
function TBSpin end

function TBSpin(Rspace::BareRspace, operator::AbstractVector)
    @assert !isempty(operator) "empty operator"
    @assert !isempty(operator[1]) "empty operator"
    @assert operator[1] isa AbstractMatrix "operator must be a matrix"
    v = operator[1][1, 1]
    @assert v isa AbstractVector && length(v) == 3 "each element must be 3-vector"
    T = real(eltype(v))
    M = Matrix{MVec3{Complex{T}}}
    return TBOperator{M}("Spin", Rspace, operator)
end

"""
    $(TYPEDEF)

A struct for interpolating tight-binding spin operator on given kpoints.

# Fields
$(FIELDS)
"""
struct SpinInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian.
    Since we interpolate on kpoints in Bloch gauge, we need to store the Hamiltonain.
    """
    hamiltonian::TBOperator

    """R-space spin operator."""
    spin::TBOperator
end

"""Interpolate the spin operator and transform it to Bloch gauge."""
function (interp::SpinInterpolator)(kpoints::AbstractVector{<:AbstractVector}; kwargs...)
    # R-space Hamiltonain
    H_R = interp.hamiltonian
    # k-space Hamiltonian
    H_k = invfourier(H_R, kpoints)
    # diagonalize
    _, gauges = eigen(H_k)

    # Wannier-gauge k-space spin operator
    Sᵂ_k = invfourier(interp.position, kpoints)
    # transform to Bloch gauge
    S_k = map(zip(Sᵂ_k, gauges)) do (Sᵂ, U)
        U' * Sᵂ * U
    end
    return S_k
end
