export read_w90_tb

"""
    read_w90_tb(seedname; force_mdrs=false)

Read `seedname_tb.dat` and `seedname_wsvec.dat`, return [`TBHamiltonian`](@ref) and [`TBPosition`](@ref) operators.

# Arguments
- `seedname`: the seedname of `seedname_tb.dat` and `seedname_wsvec.dat`

# Return
- A [`RVectors`](@ref).
- A [`TBHamiltonian`](@ref).
- A `rx` [`TBPosition`](@ref).
- A `ry` [`TBPosition`](@ref).
- A `rz` [`TBPosition`](@ref).

!!! note

    Usually, if you only need to interpolate operators (e.g., the band structure),
    that is only inverse Fourier transform [`invfourier`](@ref) is needed,
    then you don't need to provide `kpoints`.
    If in some cases, you want to generate Wannier gauge operators, then passing
    the `kpoints` allows you to construct a complete [`KRVectors`](@ref),
    and you can subsequently call [`fourier`](@ref) on Bloch gauge operators
    to get Wannier gauge operators.

!!! warning

    Since no atomic positions and labels are provided, the auto-generated
    `Brillouin.KPath` is probably wrong. It is recommended that you pass the
    `atom_positions` and `atom_labels` arguments so that this function can
    auto generate the correct `Brillouin.KPath`.
"""
function read_w90_tb(seedname::AbstractString; force_mdrs=false)
    wsvec = read_w90_wsvec(seedname * "_wsvec.dat")
    mdrs = wsvec.mdrs
    mdrs |= force_mdrs
    Rvectors = wsvec.Rvectors
    if mdrs
        Tvectors = wsvec.Tvectors
        Tdegens = wsvec.Tdegens
    end

    tbdat = read_w90_tbdat(seedname * "_tb.dat")

    lattice = tbdat.lattice

    Rvectors == tbdat.Rvectors || @error "R vecs in tb.dat and wsvec.dat are not identical"
    Rdegens = tbdat.Rdegens

    # grid is still unknown
    grid = Vec3(-1, -1, -1)  # no grid in wsvec.dat
    Rvectors_ws = RVectors(lattice, grid, Rvectors, Rdegens)

    function generate_TBBlocks(O)
        map(enumerate(O)) do (iR, o)
            rcryst = Rvectors[iR]
            rcart = lattice * rcryst
            return TBBlock(rcryst, rcart, o, o)
        end
    end

    if !mdrs
        return (;
            Rvectors=Rvectors_ws,
            H=generate_TBBlocks(tbdat.H),
            r_x=generate_TBBlocks(tbdat.r_x),
            r_y=generate_TBBlocks(tbdat.r_y),
            r_z=generate_TBBlocks(tbdat.r_z),
        )
    end

    Rvectors_mdrs = RVectorsMDRS(Rvectors_ws, Tvectors, Tdegens)

    return (;
        Rvectors=Rvectors_mdrs,
        H=mdrs_v1tov2(tbdat.H, Rvectors_mdrs),
        r_x=mdrs_v1tov2(tbdat.r_x, Rvectors_mdrs),
        r_y=mdrs_v1tov2(tbdat.r_y, Rvectors_mdrs),
        r_z=mdrs_v1tov2(tbdat.r_z, Rvectors_mdrs),
    )
end
