"""
    $(SIGNATURES)

Convert an integer to a subscript, e.g. 1 -> ₁, 2 -> ₂, etc.
"""
function to_subscript(i::Integer)
    @assert 0 <= i <= 9
    return Char(0x2080 + i)
end

function show_lattice(io, lattice::Mat3)
    @printf(io, "lattice:    [ax, ay, az] (Å)\n")
    for (i, vec) in enumerate(eachcol(lattice))
        @printf(io, "  a%s = %9.6f %9.6f %9.6f\n", to_subscript(i), vec...)
    end
end

function show_recip_lattice(io, recip_lattice::Mat3)
    @printf(io, "recip_lattice:  [bx, by, bz] (Å⁻¹)\n")
    for (i, vec) in enumerate(eachcol(recip_lattice))
        @printf(io, "  b%s = %9.6f %9.6f %9.6f\n", to_subscript(i), vec...)
    end
end
