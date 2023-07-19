using FastLapackInterface: HermitianEigenWs, decompose!
using ProgressMeter: Progress, next!

export HamiltonianRspace, HamiltonianKspace, eigen, interpolate

"""
    $(TYPEDEF)

A struct representing tight-binding Hamiltonian in R-space.

# Fields
$(FIELDS)
"""
struct HamiltonianRspace{M<:AbstractMatrix} <: AbstractOperatorRspace
    """the R-space domain (or called R-vectors) on which the operator is defined"""
    domain::BareRspaceDomain

    """The tight-binding operator defined on Rspace domain."""
    operator::Vector{M}
end

"""
    $(TYPEDEF)

A struct representing tight-binding Hamiltonian on a kpoint grid.

# Fields
$(FIELDS)
"""
struct HamiltonianKspace{K<:AbstractKpointContainer,M<:AbstractMatrix} <:
       AbstractOperatorKspace{K}
    """a [`KpointGrid`](@ref) or [`KpointList`](@ref) on which the operator is defined"""
    domain::K

    """the tight-binding operator defined on `domain`"""
    operator::Vector{M}
end

"""
    $(SIGNATURES)

Construct a `HamiltonianKspace` from `domain` and `eigenvalues`.

Construct a `HamiltonianKspace` from a Wannierization [`Model`](@ref).
"""
function HamiltonianKspace(
    domain::AbstractKpointContainer, eigenvalues::AbstractVector{T}
) where {T<:AbstractVector}
    Hᵏ = [Diagonal(εk) for εk in eigenvalues]
    return HamiltonianKspace(domain, Hᵏ)
end

function HamiltonianKspace(model::Model, gauges::AbstractVector=model.gauges)
    Hᵏ = transform_gauge(model.eigenvalues, gauges)
    return HamiltonianKspace(model.kgrid, Hᵏ)
end

"""
    $(SIGNATURES)

Construct a [`HamiltonianRspace`](@ref) from a Wannierization [`Model`](@ref).

# Arguments
- `model`: the Wannierization [`Model`](@ref)

# Keyword Arguments
- `MDRS`: whether to use MDRS interpolation
"""
function HamiltonianRspace(domain::AbstractRspaceDomain, hamiltonian_k::HamiltonianKspace)
    Hᴿ = fourier(hamiltonian_k, domain)
    bare_Rdomain, bare_Hᴿ = simplify(domain, Hᴿ)
    return HamiltonianRspace(bare_Rdomain, bare_Hᴿ)
end

function HamiltonianRspace(domain::AbstractRspaceDomain, model::Model)
    Hᵏ = HamiltonianKspace(model)
    return HamiltonianRspace(domain, Hᵏ)
end

function HamiltonianRspace(model::Model; MDRS::Bool=true)
    domain = generate_Rspace_domain(model, MDRS)
    return HamiltonianRspace(domain, model)
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

@inline function LinearAlgebra.eigen(hamiltonian::HamiltonianKspace)
    @assert length(hamiltonian) > 0 "empty hamiltonian"
    T = eltype(hamiltonian[1])
    nkpts = n_kpoints(hamiltonian)
    nwann = n_wannier(hamiltonian)

    caches = [HermitianEigenWs(zeros(T, nwann, nwann)) for _ in 1:Threads.nthreads()]
    progress = Progress(nkpts, 1, "Diagonalizing H(k)...")

    eigenvals = [zeros(real(T), nwann) for _ in 1:n_kpoints(hamiltonian)]
    eigenvecs = [zeros(T, nwann, nwann) for _ in 1:n_kpoints(hamiltonian)]

    Threads.@threads for ik in 1:nkpts
        tid = Threads.threadid()
        eigenvecs[ik] .= hamiltonian[ik]
        eigen!(eigenvals[ik], eigenvecs[ik], caches[tid])
        # this is slow
        # e = eigen(Hermitian(eigenvecs[ik]))
        # eigenvals[ik] .= e.values
        # eigenvecs[ik] .= e.vectors
        next!(progress)
    end
    return eigenvals, eigenvecs
end

"""Transform the gauge of Hamiltonian to Bloch gauge."""
function transform_gauge!(hamiltonian_k::HamiltonianKspace, gauges::AbstractVector)
    map(enumerate(zip(hamiltonian_k.operator, gauges))) do (ik, (H, U))
        hamiltonian_k.operator[ik] .= U' * H * U
    end
    return hamiltonian_k
end

function transform_gauge(hamiltonian_k::HamiltonianKspace, gauges::AbstractVector)
    UtHU = map(zip(hamiltonian_k.operator, gauges)) do (H, U)
        U' * H * U
    end
    return HamiltonianKspace(hamiltonian_k.domain, UtHU)
end

"""Interpolate the Hamiltonian operator and transform it to Bloch gauge."""
function interpolate(hamiltonian_R::HamiltonianRspace, kdomain::AbstractKpointContainer)
    Hᵏ = invfourier(hamiltonian_R, kdomain)
    hamiltonian_k = HamiltonianKspace(kdomain, Hᵏ)
    eigenvals, _ = eigen(hamiltonian_k)
    return eigenvals
end

function interpolate(hamiltonian_R::HamiltonianRspace, kpoints::AbstractVector)
    klist = KpointList(reciprocal_lattice(real_lattice(hamiltonian_R)), kpoints)
    return interpolate(hamiltonian_R, klist)
end

function interpolate(hamiltonian_R::HamiltonianRspace, kpi::KPathInterpolant)
    kpoints = get_kpoints(kpi)
    return interpolate(hamiltonian_R, kpoints)
end
