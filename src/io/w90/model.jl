export read_w90, read_w90_interp, write_w90

"""
    read_w90(seedname::AbstractString; amn=true, orthonorm_amn=true, mmn=true, eig=true)

Read `win`, and optionally `amn`, `mmn`, `eig`.

# Keyword arguments
- amn: if `true`, read `amn` file
- orthonorm_amn: Lowdin orthonormalization after reading `amn` matrices.
    Should be `true` for most cases, since usually the input `amn` matrices are
    projections onto atomic orbitals, and are not unitary or semi-unitary.
- mmn: if `true`, read `mmn` file
- eig: if `true`, read `eig` file
"""
function read_w90(
    seedname::AbstractString;
    amn::Bool=true,
    orthonorm_amn::Bool=true,
    mmn::Bool=true,
    eig::Bool=true,
)
    win = read_win("$seedname.win")

    n_bands = win.num_bands
    n_wann = win.num_wann
    kpoints = win.kpoints
    kgrid = Vec3{Int}(win.mp_grid)
    n_kpts = prod(kgrid)

    lattice = win.unit_cell_cart
    recip_lattice = get_recip_lattice(lattice)

    bvectors = get_bvectors(kpoints, recip_lattice; kmesh_tol=get(win, :kmesh_tol, 1e-6))

    n_bvecs = bvectors.n_bvecs
    if mmn
        M, kpb_k_mmn, kpb_b_mmn = read_mmn("$seedname.mmn")

        # check consistency for mmn
        (n_bands, n_bands) != size(M[1][1]) && error("n_bands != size(M[1][1])")
        n_bvecs != length(M[1]) && error("n_bvecs != length(M[1])")
        n_kpts != length(M) && error("n_kpts != length(M)")
        bvectors.kpb_k != kpb_k_mmn && error("kpb_k != kpb_k from mmn file")
        bvectors.kpb_b != kpb_b_mmn && error("kpb_b != kpb_b from mmn file")
    else
        M = [[zeros(ComplexF64, n_bands, n_bands) for ib in 1:n_bvecs] for ik in 1:n_kpts]
    end

    if amn
        if orthonorm_amn
            U = read_orthonorm_amn("$seedname.amn")
        else
            U = read_amn("$seedname.amn")
        end
        n_bands != size(U[1], 1) && error("n_bands != size(U[1], 1)")
        n_wann != size(U[1], 2) && error("n_wann != size(U[1], 2)")
        n_kpts != length(U) && error("n_kpts != length(U)")
    else
        U = [zeros(ComplexF64, n_bands, n_wann) for i in 1:n_kpts]
    end

    if eig
        E = read_eig("$seedname.eig")
        n_bands != length(E[1]) && error("n_bands != length(E[1])")
        n_kpts != length(E) && error("n_kpts != length(E)")
    else
        E = [zeros(Float64, n_bands) for i in 1:n_kpts]
    end

    if eig && n_bands != n_wann
        dis_froz_max = get(win, :dis_froz_max, nothing)
        dis_froz_min = get(win, :dis_froz_min, -Inf)
        if isnothing(dis_froz_max)
            frozen_bands = [falses(n_bands) for i in 1:n_kpts]
        else
            frozen_bands = get_frozen_bands(E, dis_froz_max, dis_froz_min)
        end
    else
        frozen_bands = [falses(n_bands) for i in 1:n_kpts]
    end

    atom_labels = map(x -> string(x.first), win.atoms_frac)
    atom_positions = map(x -> x.second, win.atoms_frac)

    return Model(
        lattice,
        atom_positions,
        atom_labels,
        kgrid,
        kpoints,
        bvectors,
        frozen_bands,
        M,
        U,
        E,
    )
end

"""
    read_w90_interp(seedname::AbstractString; chk=true, amn=nothing, mdrs=nothing)

Return an `TBHamiltonian` for Wannier interpolation.

# Keyword arguments
- chk: if `true`, read `chk` or `chk.fmt` file to get the unitary matrices,
    otherwise read `amn` file for unitary matrices.
- amn: filename for `amn`, read this `amn` file for unitary matrices.
    Only used if not reading `chk` and `amn` is given.
- mdrs: use MDRS interpolation, else Wigner-Seitz interpolation.
    If is `nothing`, detect from `win` file;
    and if no `use_ws_distance` in `win` file, default to `true`.
"""
function read_w90_interp(
    seedname::AbstractString;
    chk::Bool=true,
    amn::Union{Nothing,AbstractString}=nothing,
    mdrs::Union{Nothing,Bool}=nothing,
    kwargs...,
)
    # read for kpoint_path, use_ws_distance
    win = read_win("$seedname.win")
    if isnothing(mdrs)
        mdrs = get(win, :use_ws_distance, true)
    end

    if chk
        model = read_w90(seedname; amn=false)
        if isfile("$seedname.chk.fmt")
            fchk = read_chk("$seedname.chk.fmt")
        else
            fchk = read_chk("$seedname.chk")
        end
        model.U .= get_U(fchk)
        centers = fchk.r
    else
        if isnothing(amn)
            model = read_w90(seedname)
        else
            model = read_w90(seedname; amn=false)
            model.U .= read_orthonorm_amn(amn)
        end
        centers = center(model)
    end
    # from cartesian to fractional
    centers = map(c -> inv(model.lattice) * c, centers)

    if mdrs
        Rvecs = get_Rvectors_mdrs(model.lattice, model.kgrid, centers)
    else
        Rvecs = get_Rvectors_ws(model.lattice, model.kgrid)
    end

    hami = TBHamiltonian(model, Rvecs; kwargs...)

    if haskey(win, :kpoint_path)
        kpath = get_kpath(win.unit_cell_cart, win.kpoint_path)
    else
        # an empty kpath
        points = Dict{Symbol,Vec3{Float64}}()
        paths = Vector{Vector{Symbol}}()
        basis = ReciprocalBasis([v for v in eachcol(model.recip_lattice)])
        setting = Ref(Brillouin.LATTICE)
        kpath = KPath{3}(points, paths, basis, setting)
    end

    return hami, kpath
end

"""
    write_w90(seedname::AbstractString, model::Model)

Write `Model` into `eig`, `mmn`, `amn` files.

# Keyword arguments
- `binary`: write `eig`, `mmn`, and `amn` in Fortran binary format
"""
function write_w90(seedname::AbstractString, model::Model; binary::Bool=false)
    # Note I need to use basename! If seedname is absolute path the joinpath
    # will just return seedname.
    outname(suffix::String) = "$seedname.$suffix"

    write_eig(outname("eig"), model.E; binary=binary)

    kpb_k = model.bvectors.kpb_k
    kpb_b = model.bvectors.kpb_b
    write_mmn(outname("mmn"), model.M, kpb_k, kpb_b; binary=binary)

    write_amn(outname("amn"), model.U; binary=binary)

    return nothing
end
