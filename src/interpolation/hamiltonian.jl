using FastLapackInterface: HermitianEigenWs, decompose!
using Base.Iterators: partition
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
        eigenvalues::AbstractMatrix,
        gauges::AbstractArray{<:Complex, 3},
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
function TBHamiltonian(model::Model, gauges::AbstractArray{<:Complex, 3} = model.gauges; kwargs...)
    Rspace = generate_Rspace(model; kwargs...)
    return TBHamiltonian(Rspace, kpoints(model), model.eigenvalues, gauges)
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

function (interp::HamiltonianInterpolator)(kpoints::AbstractVector{<:AbstractVector})
    Hᵏ = invfourier(interp.hamiltonian, kpoints)
    return eigen(Hᵏ)
end

@inline function LinearAlgebra.eigen(A::AbstractMatrix, ws::HermitianEigenWs)
    return Eigen(decompose!(ws, 'V', 'A', 'U', A, 0.0, 0.0, 0, 0, 1.0e-16)...)
end

@inline function LinearAlgebra.eigen!(
        eigenvals::AbstractVector, eigenvecs::AbstractMatrix, ws::HermitianEigenWs
    )
    decompose!(ws, 'V', 'A', 'U', eigenvecs, 0.0, 0.0, 0, 0, 1.0e-16)
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
    # Partition kpoints into chunks for each tasks
    # https://julialang.org/blog/2023/07/PSA-dont-use-threadid/#better_fix_work_directly_with_tasks
    chunk_size = max(1, nkpts ÷ n_threads)
    kpt_chunks = partition(1:nkpts, chunk_size)

    progress = Progress(
        nkpts; dt = 1, desc = "Diagonalizing matrices using $n_threads threads..."
    )

    tasks = map(kpt_chunks) do chunk
        # Each chunk gets its own spawned task that does its own local,
        # sequential work and then returns the result
        Threads.@spawn begin
            cache = HermitianEigenWs(zeros(T, nwann, nwann))
            for ik in chunk
                eigenvecs[ik] .= hamiltonian[ik]
                eigen!(eigenvals[ik], eigenvecs[ik], cache)
                # this is slower
                # e = eigen(Hermitian(eigenvecs[ik]))
                # eigenvals[ik] .= e.values
                # eigenvecs[ik] .= e.vectors
                next!(progress)
            end
        end
    end
    wait.(tasks)

    return nothing
end

function LinearAlgebra.eigen!(
        eigenvals::AbstractMatrix,
        eigenvecs::AbstractArray{<:Any, 3},
        hamiltonian::AbstractVector{<:AbstractMatrix},
    )
    nkpts = length(hamiltonian)
    eigenvals_v = [view(eigenvals, :, ik) for ik in 1:nkpts]
    eigenvecs_v = [view(eigenvecs, :, :, ik) for ik in 1:nkpts]
    return eigen!(eigenvals_v, eigenvecs_v, hamiltonian)
end

function LinearAlgebra.eigen!(
        eigenvals::AbstractMatrix,
        eigenvecs::AbstractArray{<:Any, 3},
        hamiltonian::AbstractArray{<:Any, 3},
    )
    nkpts = size(hamiltonian, 3)
    ham_v = [view(hamiltonian, :, :, ik) for ik in 1:nkpts]
    return eigen!(eigenvals, eigenvecs, ham_v)
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
