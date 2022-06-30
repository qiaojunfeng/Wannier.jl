import LinearAlgebra as LA


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
    for ik = 1:n_kpts
        Uₖ = U[:, :, ik]
        # Hamiltonian in WF basis
        Hₖ = LA.Hermitian(Uₖ' * LA.diagm(0 => E[:, ik]) * Uₖ)

        # Diagonalize
        Dₖ, Vₖ = LA.eigen(Hₖ)
        Ev[:, ik] = Dₖ[1:n_val]
        Ec[:, ik] = Dₖ[n_val+1:end]

        # Although the diagonalization destroy the smoothness of gauge,
        # one can run max localization to smooth these random gauge, since
        # there is no disentanglement. However, this requires many iterations
        # of max localization.
        # I store the gauge rotation due to diagonalization.
        Vv[:, :, ik] = Vₖ[:, 1:n_val]
        Vv[:, :, ik] = Vₖ[:, n_val+1:end]
    end

    Ev, Ec, Vv, Vc
end


"""
Split MLWFs EIG, MMN, UNK files into valence and conduction groups.
"""
function rotate_mmn(M::Array{T,4}, kpb_k::Matrix{Int}, U::Array{T,3}) where {T<:Complex}
    n_bands, n_wann = size(U)
    n_kpts = size(M, 4)
    n_bvecs = size(M, 3)

    n_bands != size(M, 1) && error("incompatible n_bands")

    # Fill MMN
    N = similar(M, n_wann, n_wann, n_bvecs, n_kpts)

    for ik = 1:n_kpts
        for ib = 1:n_bvecs
            ik2 = kpb_k[ib, ik]

            U₁ = U[:, :, ik]
            U₂ = U[:, :, ik2]

            N[:, :, ib, ik] = U₁' * M[:, :, ib, ik] * U₂
        end
    end

    N
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
    outdir_val::String = "val",
    outdir_cond::String = "cond",
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
        match = match(regex, unk)
        match === nothing && continue

        ik = parse(Int, match.captures[1])

        ik2, Ψ = read_unk(unk)
        @assert ik2 == ik

        Uvₖ = Uv[:, :, ik]
        Ucₖ = Uc[:, :, ik]

        ΨUv = Ψ * Uvₖ
        ΨUc = Ψ * Ucₖ

        val = joinpath(outdir_val, unk)
        write_unk(val, ik, ΨUv)

        cond = joinpath(outdir_cond, unk)
        write_unk(cond, ik, ΨUc)

        println("ik = ", ik, " files written: ", val, " ", cond)
    end

    nothing
end


function ones_amn(T::Type, n_wann::Int, n_kpts::Int)
    A = zeros(T, n_wann, n_wann, n_kpts)
    Iₖ = LA.diagm(0 => ones(n_wann))

    for ik = 1:n_kpts
        A[:, :, ik] = Iₖ
    end

    A
end


"""
Extract AMN matrices from chk.
"""
function get_amn(chk::Chk)
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann

    U = similar(chk.U, n_bands, n_wann, n_kpts)

    if !chk.have_disentangled
        U .= chk.U
        return U
    end

    for ik = 1:n_kpts
        # Uᵈ: semi-unitary matrices from disentanglement
        # U: unitary matrices from maximal localization
        U[:, :, ik] = chk.Uᵈ[:, :, ik] * chk.U[:, :, ik]
    end

    U
end


"""
Write splitted AMN/MMN/EIG/UNK(optional) files into valence and conduction groups.

Args:
    seedname: _description_
    n_val: number of valence WFs
    return_rotation: Return rotation eigenvectors or not, useful for further
    rotation of UNK files or other operators.
"""
function split_model(model::Model, n_val::Int, return_rotation::Bool = false)
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
    for ik = 1:n_kpts
        UVv[:, :, ik] = U[:, :, ik] * Vv[:, :, ik]
        UVc[:, :, ik] = U[:, :, ik] * Vc[:, :, ik]
    end

    # MMN
    Mv = rotate_mmn(M, kpb_k, UVv)
    Mc = rotate_mmn(M, kpb_k, UVc)

    # AMN
    Av = ones_amn(eltype(M), n_val, n_kpts)
    Ac = ones_amn(eltype(M), n_wann - n_val, n_kpts)

    model_v = Model(
        model.lattice,
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
        model.kgrid,
        model.kpoints,
        model.bvectors,
        zeros(Bool, 0, 0),
        Mc,
        Ac,
        Ec,
    )

    if return_rotation
        return model_v, model_c, UVv, UVc
    end

    model_v, model_c
end


"""
Write splitted AMN/MMN/EIG/UNK(optional) files into valence and conduction groups.

Args:
    seedname: _description_
    n_val: number of valence WFs
    unk: Write UNK files. Defaults to False.
    outdir_val: _description_. Defaults to "val".
    outdir_cond: _description_. Defaults to 'cond'.
"""
function split_wannierize(model::Model, n_val::Int, return_rotation::Bool = false)
    splitted = split_model(model, n_val, return_rotation)

    if return_rotation
        model_v, model_c, UVv, UVc = splitted
    else
        model_v, model_c = splitted
    end

    model_v.A = parallel_transport(model_v)
    model_c.A = parallel_transport(model_c)

    if return_rotation
        return model_v, model_c, UVv, UVc
    end

    model_v, model_c
end
