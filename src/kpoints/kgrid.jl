export KpointGrid

"""
    $(TYPEDEF)

A uniformlly-spaced kpoint grid.

# Fields
$(FIELDS)
"""
struct KpointGrid{T<:Real} <: AbstractKpointContainer
    """reciprocal cell, 3 * 3, each column is a reciprocal lattice vector
    in Å⁻¹ unit"""
    recip_lattice::Mat3{T}

    """number of kpoints along three reciprocal lattice vectors"""
    kgrid_size::Vec3{Int}

    """kpoint fractional coordinates, length-`n_kpoints` vector of `Vec3`"""
    kpoints::Vector{Vec3{T}}
end

Base.size(kgrid::KpointGrid) = Tuple(kgrid.kgrid_size)

function KpointGrid(
    recip_lattice::AbstractMatrix, kgrid_size::AbstractVector, kpoints::AbstractVector
)
    @assert length(kpoints) > 0 "kpoints is empty"
    T = promote_type(eltype(recip_lattice), eltype(kpoints[1]))
    return KpointGrid{T}(Mat3(recip_lattice), Vec3(kgrid_size), Vec3.(kpoints))
end

function KpointGrid(recip_lattice::AbstractMatrix, kgrid_size::AbstractVector)
    @assert all(kgrid_size .> 0) "kgrid_size must be positive"
    return KpointGrid(recip_lattice, kgrid_size, get_kpoints(kgrid_size))
end

function Base.show(io::IO, ::MIME"text/plain", kgrid::KpointGrid)
    @printf(io, "kpoint type         :  %s\n\n", nameof(typeof(kgrid)))
    _print_recip_lattice(io, kgrid.recip_lattice)
    println(io)
    @printf(io, "kgrid_size  =  %d %d %d\n", size(kgrid)...)
    @printf(io, "n_kpoints   =  %d", n_kpoints(kgrid))
end
