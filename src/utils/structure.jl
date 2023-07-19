using PeriodicTable: elements
using WannierIO: get_lattice, get_recip_lattice
export real_lattice, reciprocal_lattice

"""
    $(SIGNATURES)

Get atom number from symbol.
"""
function get_atom_number(symbol::AbstractString)
    return get_atom_number([symbol])[1]
end

"""
    $(SIGNATURES)

Get atom number from symbol.
"""
function get_atom_number(symbol::AbstractVector)
    table = [e.symbol for e in elements]
    return [findfirst(x -> x == s, table) for s in symbol]
end

function real_lattice(recip_lattice::AbstractMatrix)
    return get_lattice(Mat3(recip_lattice))
end

function reciprocal_lattice(lattice::AbstractMatrix)
    return get_recip_lattice(Mat3(lattice))
end
