export read_w90_tb

"""
    read_w90_tb(seedname;els)

Read `seedname_tb.dat` and `seedname_wsvec.dat`, return [`TBHamiltonian`](@ref) and [`TBPosition`](@ref) operators.

# Arguments
- `seedname`: the seedname of `seedname_tb.dat` and `seedname_wsvec.dat`

# Return
- A [`RVectors`](@ref).
- A [`TBHamiltonian`](@ref).
- A `Rx` [`TBPosition`](@ref).
- A `Ry` [`TBPosition`](@ref).
- A `Rz` [`TBPosition`](@ref).

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
function read_w90_tb(
    seedname::AbstractString;
    force_mdrs = false
)

    mdrs, wsvec = read_w90_wsvec(seedname * "_wsvec.dat")
    mdrs |= force_mdrs 
    R = wsvec.R
    if mdrs
        T = wsvec.T
        Nᵀ = wsvec.Nᵀ
    end

    tbdat = read_w90_tbdat(seedname * "_tb.dat")
    
    lattice = tbdat.lattice
    
    R == tbdat.R || @error "R vecs in tb.dat and wsvec.dat are not identical"
    N = tbdat.N

    # grid is still unknown
    grid = Vec3(-1, -1, -1)  # no grid in wsvec.dat
    Rvecs = RVectors(lattice, grid, R, N)

    function generate_TBBlocks(O)
        map(enumerate(O)) do (iR, o)
            rcryst = R[iR]
            rcart  = lattice * rcryst
            return TBBlock(rcryst, rcart, o, o)
        end
    end
    
    if !mdrs
        return (R  = Rvecs,
                H  = generate_TBBlocks(tbdat.H),
                Rx = generate_TBBlocks(tbdat.Rx),
                Ry = generate_TBBlocks(tbdat.Ry),
                Rz = generate_TBBlocks(tbdat.Rz))
    end

    Rvecs_mdrs = RVectorsMDRS(Rvecs, T, Nᵀ)

    return (R = Rvecs_mdrs,
            H = mdrs_v1tov2(tbdat.H, Rvecs_mdrs),
            Rx = mdrs_v1tov2(tbdat.Rx, Rvecs_mdrs),
            Ry = mdrs_v1tov2(tbdat.Ry, Rvecs_mdrs),
            Rz = mdrs_v1tov2(tbdat.Rz, Rvecs_mdrs))
end
