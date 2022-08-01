using PeriodicTable: elements

"""Get atom number from symbol"""
function get_atom_number(symbol::String)
    return get_atom_number([symbol])[1]
end

"""Get atom number from symbol"""
function get_atom_number(symbol::AbstractVector{String})
    table = [e.symbol for e in elements]
    return [findfirst(x -> x == s, table) for s in symbol]
end
