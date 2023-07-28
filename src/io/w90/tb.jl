export read_w90_tb

"""
    $(SIGNATURES)

Read `prefix_tb.dat` and `prefix_wsvec.dat` and construct tight-binding models.

# Arguments
- `prefix`: the prefix of `prefix_tb.dat` and `prefix_wsvec.dat`

# Return
- a [`HamiltonianRspace`](@ref)
- a [`PositionRspace`](@ref)

!!! note
    This will call [`simplify`](@ref) to absorb the R-vector degeneracies and
    T-vectors into the operator, leading to faster interpolations.
"""
function read_w90_tb(prefix::AbstractString)
    wsvec = read_w90_wsvec(prefix * "_wsvec.dat")
    tbdat = read_w90_tbdat(prefix * "_tb.dat")
    @assert wsvec.Rvectors == tbdat.Rvectors "R-vectors in tb.dat and wsvec.dat are not identical"

    if wsvec.mdrs
        Rspace = MDRSRspace(
            tbdat.lattice, tbdat.Rvectors, tbdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens
        )
    else
        Rspace = WignerSeitzRspace(tbdat.lattice, tbdat.Rvectors, tbdat.Rdegens)
    end

    bare_Rspace, bare_H = simplify(Rspace, tbdat.H)
    hamiltonian = TBHamiltonian(bare_Rspace, bare_H)

    # convert to matrix of MVec3, here mutable since we might need to invoke
    # some in-place functions in later interpolation steps
    pos_vecs = map(zip(tbdat.r_x, tbdat.r_y, tbdat.r_z)) do (x, y, z)
        MVec3.(x, y, z)
    end
    _, bare_pos = simplify(Rspace, pos_vecs)
    position = TBPosition(bare_Rspace, bare_pos)

    return (; hamiltonian, position)
end
