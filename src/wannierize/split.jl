using LinearAlgebra

"""
Split eigenvalues into two groups.

E: eigenvalues
U: (semi-)Unitary matrices gauge transformation from Wannierization

E.g., groups for valence bands and conduction bands.
"""
function split_eig(E::Matrix{T}, U::Array{Complex{T},3}, n_val::Int) where {T<:Real}
    n_bands, n_wann, n_kpts = size(U)
    size(E, 1) != n_bands && error("incompatible n_bands")
    size(E, 2) != n_kpts && error("incompatible n_kpts")

    n_val < 1 && error("n_val < 0")
    n_val >= n_wann && error("n_val >= n_wann")

    Ev = similar(E, n_val, n_kpts)
    Ec = similar(E, n_wann - n_val, n_kpts)

    # Eigenvectors from diagonalization of Wannier Hamiltonian
    Vv = similar(U, n_wann, n_val, n_kpts)
    Vc = similar(U, n_wann, n_wann - n_val, n_kpts)

    # Since valence and conduction are splitd by a gap,
    # by diagonalizing the Wannier Hamiltonian and using the eigenvalues,
    # we can demix the WFs into two groups for valence and conduction, respectively.
    for ik in 1:n_kpts
        Uₖ = U[:, :, ik]
        # Hamiltonian in WF basis
        Hₖ = Hermitian(Uₖ' * diagm(0 => E[:, ik]) * Uₖ)

        # Diagonalize
        Dₖ, Vₖ = eigen(Hₖ)
        Ev[:, ik] = Dₖ[1:n_val]
        Ec[:, ik] = Dₖ[(n_val + 1):end]

        # Although the diagonalization destroy the smoothness of gauge,
        # one can run max localization to smooth these random gauge, since
        # there is no disentanglement. However, this requires many iterations
        # of max localization.
        # I store the gauge rotation due to diagonalization.
        Vv[:, :, ik] = Vₖ[:, 1:n_val]
        Vc[:, :, ik] = Vₖ[:, (n_val + 1):end]
    end

    return Ev, Ec, Vv, Vc
end

"""
Rotate UNK files.

These are large matrices, write to disk.

dir: directory where UNK files are stored.
Uv: the Wannier (semi-)unitary matrix for rotating valence bands.
Uc: the Wannier (semi-)unitary matrix for rotating conduction bands.
"""
function split_unk(
    dir::String,
    Uv::Array{T,3},
    Uc::Array{T,3},
    outdir_val::String="val",
    outdir_cond::String="cond",
) where {T<:Complex}
    n_kpts = size(Uv, 3)

    n_kpts != size(Uc, 3) && error("incompatible n_kpts")

    !isdir(outdir_val) && mkdir(outdir_val)
    !isdir(outdir_cond) && mkdir(outdir_cond)

    println("UNK files will be written in: ")
    println("    valence   : ", outdir_val)
    println("    conduction: ", outdir_cond)

    regex = r"UNK(\d{5})\.\d"

    for unk in readdir(dir)
        m = match(regex, unk)
        m === nothing && continue

        ik = parse(Int, m.captures[1])

        ik2, Ψ = read_unk(joinpath(dir, unk))
        @assert ik2 == ik

        Uvₖ = Uv[:, :, ik]
        Ucₖ = Uc[:, :, ik]

        # * does not support multiplication of high-dimensional arrays,
        # reshape it to matrix
        n_gx, n_gy, n_gz, n_bands = size(Ψ)
        Ψ = reshape(Ψ, :, n_bands)
        # rotate
        ΨUv = Ψ * Uvₖ
        ΨUc = Ψ * Ucₖ
        # reshape back
        ΨUv = reshape(ΨUv, n_gx, n_gy, n_gz, size(Uvₖ, 2))
        ΨUc = reshape(ΨUc, n_gx, n_gy, n_gz, size(Ucₖ, 2))

        val = joinpath(outdir_val, unk)
        write_unk(val, ik, ΨUv)

        cond = joinpath(outdir_cond, unk)
        write_unk(cond, ik, ΨUc)

        println("ik = ", ik, " files written: ", val, " ", cond)
    end

    return nothing
end

"""
Write splitted AMN/MMN/EIG/UNK(optional) files into valence and conduction groups.

Args:
    seedname: _description_
    n_val: number of valence WFs
    Rotation eigenvectors are also returned, useful for further
        rotation of UNK files or other operators.
"""
function split_model(model::Model, n_val::Int)
    n_wann = model.n_wann
    n_kpts = model.n_kpts
    n_bands = model.n_bands

    !(0 < n_val < n_wann) && error("n_val <= 0 or n_val >= n_wann")

    E = model.E
    M = model.M
    kpb_k = model.bvectors.kpb_k

    # EIG
    U = model.A
    Ev, Ec, Vv, Vc = split_eig(E, U, n_val)
    UVv = similar(U, n_bands, n_val, n_kpts)
    UVc = similar(U, n_bands, n_wann - n_val, n_kpts)
    for ik in 1:n_kpts
        UVv[:, :, ik] = U[:, :, ik] * Vv[:, :, ik]
        UVc[:, :, ik] = U[:, :, ik] * Vc[:, :, ik]
    end

    # MMN
    Mv = rotate_mmn(M, kpb_k, UVv)
    Mc = rotate_mmn(M, kpb_k, UVc)

    # AMN
    Av = eyes_amn(eltype(M), n_val, n_kpts)
    Ac = eyes_amn(eltype(M), n_wann - n_val, n_kpts)

    model_v = Model(
        model.lattice,
        model.atom_positions,
        model.atom_labels,
        model.kgrid,
        model.kpoints,
        model.bvectors,
        zeros(Bool, 0, 0),
        Mv,
        Av,
        Ev,
    )
    model_c = Model(
        model.lattice,
        model.atom_positions,
        model.atom_labels,
        model.kgrid,
        model.kpoints,
        model.bvectors,
        zeros(Bool, 0, 0),
        Mc,
        Ac,
        Ec,
    )

    return model_v, model_c, UVv, UVc
end

"""
Write splitted AMN/MMN/EIG/UNK(optional) files into valence and conduction groups.

`n_val`: number of valence WFs
return: models and rotation matrices for UNK files or other operators.
"""
function split_wannierize(model::Model, n_val::Int)
    model_v, model_c, UVv, UVc = split_model(model, n_val)

    Av, _ = parallel_transport(model_v)
    model_v.A .= Av
    Ac, _ = parallel_transport(model_c)
    model_c.A .= Ac

    return model_v, model_c, UVv, UVc
end
