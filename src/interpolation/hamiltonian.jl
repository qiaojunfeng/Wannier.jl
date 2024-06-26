using FastLapackInterface: HermitianEigenWs, decompose!
using ProgressMeter: Progress, next!

export TBHamiltonian, HamiltonianInterpolator, eigen

"""Construct a tight-binding Hamiltonain in Rspace.

From a Wannierization [`Model`](@ref)."""
function TBHamiltonian end

function TBHamiltonian(Rspace::BareRspace, operator::AbstractVector)
    @assert !isempty(operator) "empty operator"
    T = real(eltype(operator[1]))
    M = Matrix{Complex{T}}
    return TBOperator{M}("Hamiltonian", Rspace, operator)
end

function TBHamiltonian(
    Rspace::AbstractRspace,
    kpoints::AbstractVector,
    eigenvalues::AbstractVector,
    gauges::AbstractVector,
)
    Hᵏ = transform_gauge(eigenvalues, gauges)
    Hᴿ = fourier(kpoints, Hᵏ, Rspace)
    bare_Rspace, bare_H = simplify(Rspace, Hᴿ)
    return TBHamiltonian(bare_Rspace, bare_H)
end

"""
    $(SIGNATURES)

Construct a [`TBHamiltonian`](@ref) from a Wannierization [`Model`](@ref).

# Arguments
- `model`: the Wannierization [`Model`](@ref)

# Keyword Arguments
- `MDRS`: whether to use MDRS interpolation
"""
function TBHamiltonian(model::Model, gauges::AbstractVector=model.gauges; kwargs...)
    Rspace = generate_Rspace(model; kwargs...)
    return TBHamiltonian(Rspace, model.kpoints, model.eigenvalues, gauges)
end

"""
    $(TYPEDEF)

A struct for interpolating tight-binding Hamiltonian on given kpoints.

# Fields
$(FIELDS)
"""
struct HamiltonianInterpolator <: AbstractTBInterpolator
    """R-space Hamiltonian"""
    hamiltonian::TBOperator
end

"""Interpolate the Hamiltonian operator and transform it to Bloch gauge."""
function (interp::HamiltonianInterpolator)(kpoints::AbstractVector{<:AbstractVector})
    # to also handle `KPathInterpolant`
    kpoints = get_kpoints(kpoints)
    Hᵏ = invfourier(interp.hamiltonian, kpoints)
    return eigen(Hᵏ)
end

@inline function LinearAlgebra.eigen(A::AbstractMatrix, ws::HermitianEigenWs)
    return Eigen(decompose!(ws, 'V', 'A', 'U', A, 0.0, 0.0, 0, 0, 1e-16)...)
end

@inline function LinearAlgebra.eigen!(
    eigenvals::AbstractVector, eigenvecs::AbstractMatrix, ws::HermitianEigenWs
)
    decompose!(ws, 'V', 'A', 'U', eigenvecs, 0.0, 0.0, 0, 0, 1e-16)
    eigenvals .= ws.w
    eigenvecs .= ws.Z
    return nothing
end

function LinearAlgebra.eigen!(
    eigenvals::AbstractVector{<:AbstractVector},
    eigenvecs::AbstractVector{<:AbstractMatrix},
    hamiltonian::AbstractVector{<:AbstractMatrix},
)
    @assert length(hamiltonian) > 0 "empty hamiltonian"
    T = eltype(hamiltonian[1])
    nkpts = length(hamiltonian)
    nwann = size(hamiltonian[1], 1)

    n_threads = Threads.nthreads()
    caches = [HermitianEigenWs(zeros(T, nwann, nwann)) for _ in 1:n_threads]
    progress = Progress(nkpts, 1, "Diagonalizing matrices using $n_threads threads...")

    Threads.@threads for ik in 1:nkpts
        tid = Threads.threadid()
        eigenvecs[ik] .= hamiltonian[ik]
        eigen!(eigenvals[ik], eigenvecs[ik], caches[tid])
        # this is slower
        # e = eigen(Hermitian(eigenvecs[ik]))
        # eigenvals[ik] .= e.values
        # eigenvecs[ik] .= e.vectors
        next!(progress)
    end
    return nothing
end

function LinearAlgebra.eigen(hamiltonian::AbstractVector{<:AbstractMatrix})
    nkpts = length(hamiltonian)
    @assert nkpts > 0 "empty hamiltonian"
    nwann = size(hamiltonian[1], 1)
    T = eltype(hamiltonian[1])
    eigenvals = zeros_eigenvalues(real(T), nkpts, nwann)
    eigenvecs = zeros_gauge(T, nkpts, nwann)
    eigen!(eigenvals, eigenvecs, hamiltonian)
    return eigenvals, eigenvecs
end
