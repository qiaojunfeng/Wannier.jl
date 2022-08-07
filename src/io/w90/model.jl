
@doc raw"""
read win, and optionally amn, mmn, eig

orthonorm_amn: Lowdin orthonormalization of AMN matrices.
    Should be true for most cases, since usually the input AMN matrices are
    projections onto atomic orbitals, and are not unitary or semi-unitary.
"""
function read_w90(
    seedname::String;
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

    lattice = win.unit_cell
    recip_lattice = get_recip_lattice(lattice)

    if ismissing(win.kmesh_tol)
        bvectors = get_bvectors(kpoints, recip_lattice)
    else
        bvectors = get_bvectors(kpoints, recip_lattice; kmesh_tol=win.kmesh_tol)
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
            A = read_orthonorm_amn("$seedname.amn")
        else
            A = read_amn("$seedname.amn")
        end
        n_bands != size(A)[1] && error("n_bands != size(A)[1]")
        n_wann != size(A)[2] && error("n_wann != size(A)[2]")
        n_kpts != size(A)[3] && error("n_kpts != size(A)[3]")
    else
        A = zeros(ComplexF64, n_bands, n_wann, n_kpts)
    end

    if eig
        E = read_eig("$seedname.eig")
        n_bands != size(E)[1] && error("n_bands != size(E)[1]")
        n_kpts != size(E)[2] && error("n_kpts != size(E)[2]")
    else
        E = zeros(Float64, n_bands, n_kpts)
    end

    if eig && n_bands != n_wann
        dis_froz_max = win.dis_froz_max
        dis_froz_min = win.dis_froz_min
        if dis_froz_max !== missing
            if dis_froz_min === missing
                dis_froz_min = -Inf
            end
            frozen_bands = get_frozen_bands(E, dis_froz_max, dis_froz_min)
        else
            frozen_bands = falses(n_bands, n_kpts)
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
        A,
        E,
    )
end

"""
Return an `InterpolationModel` for Wanier interpolation.

Need to read chk file.

Rvectors: which kind of interpolation? Wigner-Seitz, MDRS, or detect from win file.
"""
function read_w90_post(seedname::String; mdrs::Union{Nothing,Bool}=nothing)
    # read for kpoint_path, use_ws_distance
    win = read_win("$seedname.win")
    if isnothing(mdrs)
        if !ismissing(win.use_ws_distance) && win.use_ws_distance
            mdrs = true
        else
            mdrs = false
        end
    end

    model = read_w90(seedname; amn=false)
    chk = read_chk("$seedname.chk.fmt")
    model.A .= get_A(chk)
    centers = chk.r
    # from cartesian to fractional
    centers = inv(model.lattice) * centers

    if mdrs
        Rvecs = get_Rvectors_mdrs(model.lattice, model.kgrid, centers)
    else
        Rvecs = get_Rvectors_ws(model.lattice, model.kgrid)
    end
    kRvecs = KRVectors(model.lattice, model.kgrid, model.kpoints, Rvecs)

    if !ismissing(win.kpoint_path)
        kpath = win.kpoint_path
    else
        # an empty kpath
        points = Dict{Symbol,Vec3{Float64}}()
        paths = Vector{Vector{Symbol}}()
        basis = ReciprocalBasis([v for v in eachcol(model.recip_lattice)])
        setting = Ref(Brillouin.LATTICE)
        kpath = KPath{3}(points, paths, basis, setting)
    end

    return InterpolationModel(model, kRvecs, kpath)
end
