export read_w90_tb

"""
    read_w90_tb(seedname::AbstractString)

Read `seedname_tb.dat` and `seedname_wsvec.dat`.

# Return
- R vectors
- Hamiltonian
- position operator
"""
function read_w90_tb(seedname::AbstractString)
    mdrs, wsvec = read_w90_wsvec(seedname * "_wsvec.dat")
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
    grid = [-1, -1, -1]  # no grid in wsvec.dat
    Rvecs = RVectors(lattice, grid, R, N)
    if !mdrs
        return Rvecs, tbdat.H, tbdat.r
    end

    Rvecs_mdrs = RVectorsMDRS(Rvecs, T, Nᵀ)
    return Rvecs_mdrs, tbdat.H, tbdat.r
end
