export KpointList

"""
    $(TYPEDEF)

A list of kpoint coordinates, can be nonuniformly distributed.

# Fields
Using these accronyms,
    - `n_kpoints`: number of kpoints
the fields are defined as follows:
$(FIELDS)
"""
struct KpointList{T<:Real} <: AbstractKpointContainer
    """reciprocal cell, 3 * 3, each column is a reciprocal lattice vector
    in Å⁻¹ unit"""
    recip_lattice::Mat3{T}

    """kpoint fractional coordinates, length-`n_kpoints` vector of `Vec3`"""
    kpoints::Vector{Vec3{T}}
end

function KpointList(recip_lattice::AbstractMatrix, kpoints::AbstractVector)
    @assert length(kpoints) > 0 "kpoints is empty"
    T = promote_type(eltype(recip_lattice), eltype(kpoints[1]))
    return KpointList{T}(Mat3(recip_lattice), Vec3.(kpoints))
end

function Base.show(io::IO, ::MIME"text/plain", klist::KpointList)
    @printf(io, "kpoint type         :  %s\n\n", nameof(typeof(klist)))
    _print_recip_lattice(io, klist.recip_lattice)
    println(io)
    @printf(io, "n_kpoints   =  %d", n_kpoints(klist))
end
