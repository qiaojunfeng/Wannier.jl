"""
Abstract type storing a set of kpoints.

# Fields
Must implement the following fields, see [`KpointGrid`](@ref) for an example:
- `recip_lattice`: the reciprocal lattice, column-major in Å⁻¹ unit
- `kpoints`: a length-`n_kpoints` vector, each elment is a fractional coordinates
"""
abstract type AbstractKpointContainer end

n_kpoints(kcontainer::AbstractKpointContainer) = length(kcontainer.kpoints)
reciprocal_lattice(kcontainer::AbstractKpointContainer) = kcontainer.recip_lattice

function Base.isapprox(a::AbstractKpointContainer, b::AbstractKpointContainer; kwargs...)
    for f in propertynames(a)
        va = getfield(a, f)
        vb = getfield(b, f)

        if va isa Vector
            all(isapprox.(va, vb; kwargs...)) || return false
        else
            isapprox(va, vb; kwargs...) || return false
        end
    end
    return true
end

"""Index using `i` of kpoints"""
function Base.getindex(kcontainer::AbstractKpointContainer, i::Integer)
    return kcontainer.kpoints[i]
end
function Base.getindex(kcontainer::AbstractKpointContainer, r::UnitRange{<:Integer})
    return kcontainer.kpoints[r]
end

Base.lastindex(kcontainer::AbstractKpointContainer) = lastindex(kcontainer.kpoints)
Base.length(kcontainer::AbstractKpointContainer) = length(kcontainer.kpoints)

function Base.iterate(kcontainer::AbstractKpointContainer, state=1)
    if state > length(kcontainer.kpoints)
        return nothing
    else
        return (kcontainer.kpoints[state], state + 1)
    end
end
