export TBSpin, SpinInterpolator, SpinProjectionInterpolator

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

function TBSpin(
    Rspace::AbstractRspace,
    kpoints::AbstractVector,
    spins::AbstractVector,
    gauges::AbstractVector,
)
    Sᵏ = transform_gauge(spins, gauges)
    Sᴿ = fourier(kpoints, Sᵏ, Rspace)
    bare_Rspace, bare_S = simplify(Rspace, Sᴿ)
    return TBSpin(bare_Rspace, bare_S)
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
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    # R-space Hamiltonain
    H_R = interp.hamiltonian
    # k-space Hamiltonian
    H_k = invfourier(H_R, kpoints)
    # diagonalize
    _, gauges = eigen(H_k)

    # Wannier-gauge k-space spin operator
    Sᵂ_k = invfourier(interp.spin, kpoints)
    # transform to Bloch gauge
    S_k = map(zip(Sᵂ_k, gauges)) do (Sᵂ, U)
        U' * Sᵂ * U
    end
    return S_k
end

"""
    $(TYPEDEF)

A struct for interpolating tight-binding spin operator projected on a given
axis, on given kpoints.

# Fields
$(FIELDS)
"""
struct SpinProjectionInterpolator <: AbstractTBInterpolator
    spin_interpolator::SpinInterpolator

    """spin quantization axis polar angle, in radian."""
    θ::Float64

    """spin quantization axis azimuthal angle, in radian."""
    ϕ::Float64
end

n_wannier(interp::SpinProjectionInterpolator) = n_wannier(interp.spin_interpolator)
n_Rvectors(interp::SpinProjectionInterpolator) = n_Rvectors(interp.spin_interpolator)
real_lattice(interp::SpinProjectionInterpolator) = real_lattice(interp.spin_interpolator)

function SpinProjectionInterpolator(hamiltonian::TBOperator, spin::TBOperator, θ, ϕ)
    spin_interpolator = SpinInterpolator(hamiltonian, spin)
    return SpinProjectionInterpolator(spin_interpolator, θ, ϕ)
end

"""
Interpolate the spin operator and transform it to Bloch gauge.

# Keyword arguments
- `truncate`: if abs(spin) > 1, truncate them to inside [-1, 1]
"""
function (interp::SpinProjectionInterpolator)(
    kpoints::AbstractVector{<:AbstractVector}; truncate::Bool=true
)
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    S_k = interp.spin_interpolator(kpoints)

    # instead of matrix, return real part of diagonal elements
    # also truncate values > 1
    ax = Vec3(sin(interp.θ) * cos(interp.ϕ), sin(interp.θ) * sin(interp.ϕ), cos(interp.θ))
    return map(S_k) do S
        S_diag = real(diag(S))
        S_ax = map(S_diag) do svec # each band
            svec ⋅ ax
        end
        if truncate
            S_ax = clamp.(S_ax, -1, 1)
        end
        return S_ax
    end
end
