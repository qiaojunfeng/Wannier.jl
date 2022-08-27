using PeriodicTable: elements

"""
    get_atom_number(symbol::AbstractString)

Get atom number from symbol.
"""
function get_atom_number(symbol::AbstractString)
    return get_atom_number([symbol])[1]
end

"""
    get_atom_number(symbol::Vector{String})

Get atom number from symbol.
"""
function get_atom_number(symbol::AbstractVector{T}) where {T<:AbstractString}
    table = [e.symbol for e in elements]
    return [findfirst(x -> x == s, table) for s in symbol]
end
