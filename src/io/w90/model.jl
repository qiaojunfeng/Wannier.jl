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

    if haskey(win, :kmesh_tol)
        bvectors = get_bvectors(kpoints, recip_lattice; kmesh_tol=win.kmesh_tol)
    else
        bvectors = get_bvectors(kpoints, recip_lattice)
    end
    n_bvecs = bvectors.n_bvecs

    if mmn
        M, kpb_k_mmn, kpb_b_mmn = read_mmn("$seedname.mmn")

        # check consistency for mmn
        n_bands != size(M)[1] && error("n_bands != size(M)[1]")
        n_bvecs != size(M)[3] && error("n_bvecs != size(M)[3]")
        n_kpts != size(M)[4] && error("n_kpts != size(M)[4]")
        bvectors.kpb_k != kpb_k_mmn && error("kpb_k != kpb_k from mmn file")
        bvectors.kpb_b != kpb_b_mmn && error("kpb_b != kpb_b from mmn file")
    else
        M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    end

    if amn
        if orthonorm_amn
            U = read_orthonorm_amn("$seedname.amn")
        else
            U = read_amn("$seedname.amn")
        end
        n_bands != size(U)[1] && error("n_bands != size(U)[1]")
        n_wann != size(U)[2] && error("n_wann != size(U)[2]")
        n_kpts != size(U)[3] && error("n_kpts != size(U)[3]")
    else
        U = zeros(ComplexF64, n_bands, n_wann, n_kpts)
    end

    if eig
        E = read_eig("$seedname.eig")
        n_bands != size(E)[1] && error("n_bands != size(E)[1]")
        n_kpts != size(E)[2] && error("n_kpts != size(E)[2]")
    else
        E = zeros(Float64, n_bands, n_kpts)
    end

    if eig && n_bands != n_wann
        dis_froz_max = get(win, :dis_froz_max, nothing)
        dis_froz_min = get(win, :dis_froz_min, -Inf)
        if isnothing(dis_froz_max)
            frozen_bands = falses(n_bands, n_kpts)
        else
            frozen_bands = get_frozen_bands(E, dis_froz_max, dis_froz_min)
        end
    else
        frozen_bands = falses(n_bands, n_kpts)
    end

    return Model(
        lattice,
        win.atoms_frac,
        win.atom_labels,
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

Return an `InterpModel` for Wannier interpolation.

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
    centers = inv(model.lattice) * centers

    if mdrs
        Rvecs = get_Rvectors_mdrs(model.lattice, model.kgrid, centers)
    else
        Rvecs = get_Rvectors_ws(model.lattice, model.kgrid)
    end
    kRvecs = KRVectors(model.lattice, model.kgrid, model.kpoints, Rvecs)

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

    Hᵏ = get_Hk(model.E, model.U)
    Hᴿ = fourier(kRvecs, Hᵏ)

    return InterpModel(kRvecs, kpath, Hᴿ)
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
