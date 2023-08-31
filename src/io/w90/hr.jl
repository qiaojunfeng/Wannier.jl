export read_w90_hr, write_w90_hr

"""
    $(SIGNATURES)

Read `prefix_hr.dat` and `prefix_wsvec.dat` and construct tight-binding models.

# Arguments
- `prefix`: the prefix of `prefix_hr.dat` and `prefix_wsvec.dat`. If the latter
    does not exist, it will be ignored.
- `lattice`: the lattice vectors, column-wise in Å

# Return
- a [`TBHamiltonian`](@ref)

!!! note
    This will call [`simplify`](@ref) to absorb the R-vector degeneracies and
    T-vectors into the operator, leading to faster interpolations.
"""
function read_w90_hr(prefix::AbstractString, lattice::AbstractMatrix)
    dat = _raw_read_w90_hr(prefix, lattice)
    Rspace = dat.Rspace
    H = dat.hamiltonian

    bare_Rspace, bare_H = simplify(Rspace, H)
    hamiltonian = TBHamiltonian(bare_Rspace, bare_H)

    return hamiltonian
end

"""Only read hr files, without further processing"""
function _raw_read_w90_hr(prefix::AbstractString, lattice::AbstractMatrix)
    hrdat = read_w90_hrdat(prefix * "_hr.dat")
    if isfile(prefix * "_hr.dat")
        wsvec = read_w90_wsvec(prefix * "_wsvec.dat")
        @assert wsvec.Rvectors ≈ hrdat.Rvectors "R-vectors in hr.dat and wsvec.dat are not identical"
        if wsvec.mdrs
            Rspace = MDRSRspace(
                lattice, hrdat.Rvectors, hrdat.Rdegens, wsvec.Tvectors, wsvec.Tdegens
            )
        else
            Rspace = WignerSeitzRspace(lattice, hrdat.Rvectors, hrdat.Rdegens)
        end
    else
        Rspace = WignerSeitzRspace(lattice, hrdat.Rvectors, hrdat.Rdegens)
    end

    return (; Rspace, hamiltonian=hrdat.H)
end

"""
    $(SIGNATURES)

Write a tight-binding model of Hamiltonian into
`prefix_hr.dat` and `prefix_wsvec.dat` files.

# Keyword Arguments
- `skip_wsvec`: whether to skip writing `prefix_wsvec.dat`. Since the `hamiltonian.Rspace`
    is always `BareRspace`, the `prefix_wsvec.dat` file is not necessary.
"""
function write_w90_hr(prefix::AbstractString, hamiltonian::TBOperator; skip_wsvec=false)
    # the operators are always BareRspace
    if !skip_wsvec
        WannierIO.write_w90_wsvec(
            prefix * "_wsvec.dat";
            hamiltonian.Rspace.Rvectors,
            n_wann=n_wannier(hamiltonian),
        )
    end

    return WannierIO.write_w90_hrdat(
        prefix * "_hr.dat";
        hamiltonian.Rspace.Rvectors,
        Rdegens=ones(n_Rvectors(hamiltonian)),
        H=hamiltonian.operator,
    )
end
