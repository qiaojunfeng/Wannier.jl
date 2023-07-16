export read_w90, read_w90_with_chk, write_w90

"""
    read_w90(prefix; amn=true, orthonorm_amn=true, mmn=true, eig=true)

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
    prefix::AbstractString;
    amn::Bool=true,
    orthonorm_amn::Bool=true,
    mmn::Bool=true,
    eig::Bool=true,
)
    win = read_win("$prefix.win")

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
        M, kpb_k_mmn, kpb_G_mmn = read_mmn("$prefix.mmn")

        # check consistency for mmn
        (n_bands, n_bands) != size(M[1][1]) && error("n_bands != size(M[1][1])")
        n_bvecs != length(M[1]) && error("n_bvecs != length(M[1])")
        n_kpts != length(M) && error("n_kpts != length(M)")
        bvectors.kpb_k != kpb_k_mmn && error("kpb_k != kpb_k from mmn file")
        bvectors.kpb_G != kpb_G_mmn && error("kpb_G != kpb_G from mmn file")
    else
        M = [[zeros(ComplexF64, n_bands, n_bands) for ib in 1:n_bvecs] for ik in 1:n_kpts]
    end

    if amn
        if orthonorm_amn
            U = read_orthonorm_amn("$prefix.amn")
        else
            U = read_amn("$prefix.amn")
        end
        @assert n_kpts == length(U) "amn file has different n_kpts than win file: $(length(U)) != $n_kpts"
        @assert (n_bands, n_wann) == size(U[1]) "U[1] is not a n_bands x n_wann matrix: $(size(U[1])) != ($n_bands, $n_wann)"
    else
        U = [zeros(ComplexF64, n_bands, n_wann) for i in 1:n_kpts]
    end

    if eig
        E = read_eig("$prefix.eig")
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

@doc raw"""
    read_w90_with_chk(prefix, chk="$prefix.chk")

Return a `Model` with U matrix filled by that from a chk file.

# Arguments
- chk: path of chk file to get the unitary matrices.
"""
function read_w90_with_chk(prefix::AbstractString, chk::AbstractString="$prefix.chk")
    model = read_w90(prefix; amn=false)
    fchk = read_chk(chk)
    model.U .= get_U(fchk)
    return model
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
    kpb_G = model.bvectors.kpb_G
    write_mmn(outname("mmn"), model.M, kpb_k, kpb_G; binary=binary)

    write_amn(outname("amn"), model.U; binary=binary)

    return nothing
end
