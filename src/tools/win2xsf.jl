using LinearAlgebra

"""
    $(SIGNATURES)

Save crystal structure in W90 input `win` file to a `xsf` file.

# Arguments
- `xsf`: The output filename to write the crystal structure to.
- `win`: The W90 input file
"""
function win2xsf(xsf::AbstractString, win::Union{NamedTuple, AbstractDict})
    cell = get(win, :unit_cell_cart, nothing)
    isnothing(cell) && error("No unit cell found in the W90 input file.")

    if haskey(win, :atoms_frac)
        atoms_frac = win.atoms_frac
    elseif haskey(win, :atoms_cart)
        invcell = inv(cell)
        atoms_frac = map(win.atoms_cart) do p
            first(p) => invcell * last(p)
        end
    else
        error("No atoms found in the W90 input file.")
    end
    atom_numbers = get_atom_number.(first.(atoms_frac))
    atom_positions = last.(atoms_frac)

    return write_xsf(xsf, cell, atom_positions, atom_numbers)
end

function win2xsf(xsf::AbstractString, win::AbstractString)
    return win2xsf(xsf, read_win(win))
end
