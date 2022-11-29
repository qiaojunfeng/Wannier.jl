export read_w90_tb

"""
    _read_w90_tb(seedname::AbstractString)

Read `seedname_tb.dat` and `seedname_wsvec.dat`.

# Return
- R vectors
- Hamiltonian
- position operator
"""
function _read_w90_tb(seedname::AbstractString)
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

"""
    read_w90_tb(seedname; kpoints, atom_positions, atom_labels)

Read `seedname_tb.dat` and `seedname_wsvec.dat`, return an [`InterpModel`](@ref).

# Arguments
- `seedname`: the seedname of `seedname_tb.dat` and `seedname_wsvec.dat`

# Keyword Arguments
- `kpoints`: the kpoints used in Wannierization, each column is a fractional coordinate.
- `atom_positions`: columns are the fractional coordinates of atoms
- `atom_labels`: labels of atoms

# Return
- An [`InterpModel`](@ref).

!!! note

    Usually, if you only need to interpolate operators (e.g., the band structure),
    that is only inverse Fourier transform [`invfourier`](@ref) is needed,
    then you don't need to provide `kpoints`.
    If in some cases, you want to generate Wannier gauge operators, then passing
    the `kpoints` allows you to construct a complete [`KRVectors`](@ref),
    and you can subsequently call [`fourier`](@ref) on Bloch gauge operators
    to get Wannier gauge operators.

!!! warning

    Since no atomic positions and labels are provided, the auto-generated `KPath`
    is probably wrong.
    It is recommended that you pass the `atom_positions` and `atom_labels` arguments
    so that this function can auto generate the correct [`KPath`](@ref).
"""
function read_w90_tb(
    seedname::AbstractString;
    kpoints::Union{Nothing,AbstractMatrix}=nothing,
    atom_positions::Union{Nothing,AbstractMatrix}=nothing,
    atom_labels::Union{Nothing,AbstractVector{String}}=nothing,
)
    Rvecs, H, r = _read_w90_tb(seedname)
    if kpoints === nothing
        kRvecs = KRVectors(Rvecs)
    else
        size(kpoints, 1) == 3 || error("kpoints should be 3 * n_kpts")
        kgrid = Vec3{Int}(get_kgrid(kpoints))
        kRvecs = KRVectors(Rvecs.lattice, kgrid, kpoints, Rvecs)
    end

    if atom_positions === nothing || atom_labels === nothing
        # generate a fake atom
        atom_positions = zeros(3, 1)
        atom_labels = ["H"]
    end
    kpath = get_kpath(Rvecs.lattice, atom_positions, atom_labels)

    n_wann = size(H, 1)
    n_rvecs = size(H, 3)
    S = Array{ComplexF64,4}(undef, n_wann, n_wann, n_rvecs, 3)

    return InterpModel(kRvecs, kpath, H, r, S)
end
