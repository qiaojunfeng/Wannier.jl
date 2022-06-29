import LinearAlgebra as LA


"""
Split eigenvalues into two groups.

E.g., groups for valence bands and conduction bands.
"""
function split_eig(E::Matrix{T}, chk::Chk, n_val::Int) where {T<:Real}
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann

    size(E, 1) != n_bands && error("incompatible n_bands")
    size(E, 2) != n_kpts && error("incompatible n_kpts")
    n_val < 1 && error("n_val < 0")
    n_val >= n_wann && error("n_val >= n_wann")

    Ev = similar(E, n_val, n_kpts)
    Ec = similar(E, n_wann - n_val, n_kpts)

    # (semi-)Unitary matrices from Wannierization
    U = similar(chk.U, n_bands, n_wann, n_kpts)
    # Eigenvectors from diagonalization of Wannier Hamiltonian
    V = similar(U, n_wann, n_wann, n_kpts)

    # Since valence and conduction are splitd by a gap,
    # by diagonalizing the Wannier Hamiltonian and using the eigenvalues,
    # we can demix the WFs into two groups for valence and conduction, respectively.
    for ik = 1:n_kpts
        # Uᵈ: semi-unitary matrices from disentanglement
        # U: unitary matrices from maximal localization
        Wₖ = chk.Uᵈ[:, :, ik] * chk.U[:, :, ik]
        U[:, :, ik] = Wₖ

        # Hamiltonian in WF basis
        Hₖ = LA.Hermitian(Wₖ' * LA.diagm(0 => E[:, ik]) * Wₖ)
        # Diagonalize
        Dₖ, Vₖ = LA.eigen(Hₖ)
        Ev[:, ik] = Dₖ[1:n_val]
        Ec[:, ik] = Dₖ[n_val+1:end]

        # Although the diagonalization destroy the smoothness of gauge,
        # one can run max localization to smooth these random gauge, since
        # there is no disentanglement. However, this requires many iterations
        # of max localization.
        # I store the gauge rotation due to diagonalization.
        V[:, :, ik] = Vₖ
    end

    Ev, Ec, U, V
end


"""
Split MLWFs EIG, MMN, UNK files into valence and conduction groups.
"""
function split_mmn(
    M::Array{T,4},
    kpb_k::Matrix{Int},
    U::Array{T,3},
    V::Array{T,3},
    n_val::Int,
) where {T<:Complex}
    n_bands, n_wann = size(U)
    n_kpts = size(M, 4)
    n_bvecs = size(M, 3)

    n_bands != size(M, 1) && error("incompatible n_bands")
    n_wann != size(V, 1) && error("incompatible n_wann")
    n_val < 1 && error("num_val < 1")
    n_val >= n_wann && error("num_val >= n_wann")

    # Fill MMN
    Mv = similar(M, n_val, n_val, n_bvecs, n_kpts)
    Mc = similar(M, n_wann - n_val, n_wann - n_val, n_bvecs, n_kpts)

    for ik = 1:n_kpts
        for ib = 1:n_bvecs
            ik2 = kpb_k[ib, ik]

            U₁ = U[:, :, ik]
            U₂ = U[:, :, ik2]

            V₁v = V[:, 1:n_val, ik]
            V₂v = V[:, 1:n_val, ik2]

            # gauge overlap matrix
            MU = U₁' * M[:, :, ib, ik] * U₂
            Mv[:, :, ib, ik] = V₁v' * MU * V₂v

            V₁c = V[:, n_val+1:end, ik]
            V₂c = V[:, n_val+1:end, ik2]

            Mc[:, :, ib, ik] = V₁c' * MU * V₂c
        end
    end

    Mv, Mc
end


"""
Rotate UNK files.

These are large matrices, write to disk.

dir: directory where UNK files are stored.
U: the Wannier (semi-)unitary matrix for (disentanglement and) maximal localization.
V: additional rotation by user.
n_val: number of WFs for valence band.
"""
function split_unk(
    dir::String,
    U::Array{T,3},
    V::Array{T,3},
    n_val::Int;
    outdir_val::String = "val",
    outdir_cond::String = "cond",
) where {T<:Complex}

    n_kpts = size(U, 3)
    n_wann = size(U, 2)

    n_kpts != size(V, 3) && error("incompatible n_kpts")
    n_wann != size(V, 1) && error("incompatible n_wann")

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

        Uₖ = U[:, :, ik]
        Vₖv = V[:, 1:n_val, ik]
        Vₖc = V[:, n_val+1:end, ik]

        ΨU = Ψ * Uₖ

        ΨUVv = ΨU * Vₖv
        ΨUVc = ΨU * Vₖc

        val = joinpath(outdir_val, unk)
        write_unk(val, ik, ΨUVv)

        cond = joinpath(outdir_cond, unk)
        write_unk(cond, ik, ΨUVc)

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
split AMN/MMN/EIG/UNK(optional) files into valence and conduction groups.

Args:
    seedname: _description_
    n_val: number of valence WFs
    unk: Write UNK files. Defaults to False.
    vmn: Write rotation eigenvectors into AMN format. Defaults to False.
    outdir_val: _description_. Defaults to "val".
    outdir_cond: _description_. Defaults to 'cond'.
"""
function split_valence_conduction(
    seedname::String,
    n_val::Int;
    unk::Bool = false,
    vmn::Bool = false,
    outdir_val::String = "val",
    outdir_cond::String = "cond",
)
    !isdir(outdir_val) && mkdir(outdir_val)
    !isdir(outdir_cond) && mkdir(outdir_cond)

    chk = read_chk("$seedname.chk.fmt")
    E = read_eig("$seedname.eig")
    M, kpb_k, kpb_b = read_mmn("$seedname.mmn")

    n_wann = chk.n_wann
    n_kpts = chk.n_kpts

    # Note I need to use basename! If seedname is absolute path the joinpath
    # will just return seedname.
    out_val(suffix::String) = joinpath(outdir_val, basename("$seedname.$suffix"))
    out_cond(suffix::String) = joinpath(outdir_cond, basename("$seedname.$suffix"))

    # EIG
    Ev, Ec, U, V = split_eig(E, chk, n_val)
    write_eig(out_val("eig"), Ev)
    write_eig(out_cond("eig"), Ec)

    # VMN
    if vmn
        header = "Created by split_valence_conduction (valence)"
        Vv = V[:, 1:n_val, :]
        write_amn(out_val("vmn"), Vv, header)

        header = "Created by split_valence_conduction (conduction)"
        Vc = V[:, n_val+1:end, :]
        write_amn(out_cond("vmn"), Vc, header)
    end

    # MMN
    Mv, Mc = split_mmn(M, kpb_k, U, V, n_val)
    header = "Created by split_valence_conduction (valence)"
    write_mmn(out_val("mmn"), Mv, kpb_k, kpb_b, header)
    header = "Created by split_valence_conduction (conduction)"
    write_mmn(out_cond("mmn"), Mc, kpb_k, kpb_b, header)

    # AMN
    Av = ones_amn(eltype(M), n_val, n_kpts)
    header = "Created by split_valence_conduction (valence)"
    write_amn(out_val("amn"), Av, header)

    Ac = ones_amn(eltype(M), n_wann - n_val, n_kpts)
    header = "Created by split_valence_conduction (conduction)"
    write_amn(out_cond("amn"), Ac, header)

    # UNK
    if unk
        dir = dirname(seedname)
        split_unk(dir, U, V, n_val; outdir_val=outdir_val, outdir_cond=outdir_cond)
    end

    nothing
end


"""
Generate valence only MMN, EIG files from a val+cond NSCF calculation.

Args:
    seedname: _description_
    outdir: the folder for writing MMN, EIG files.
"""
function truncate_mmn_eig(
    seedname::String,
    keep_bands::AbstractVector{Int};
    outdir::String = "truncate",
)
    !isdir(outdir) && mkdir(outdir)

    E = read_eig("$seedname.eig")
    M, kpb_k, kpb_b = read_mmn("$seedname.mmn")

    E1 = E[keep_bands, :]
    write_eig(joinpath(outdir, "$seedname.eig"), E1)

    M1 = M[keep_bands, keep_bands, :, :]
    write_mmn(joinpath(outdir, "$seedname.mmn"), M1, kpb_k, kpb_b)

    nothing
end


"""
Truncate UNK files for specified bands.

Args:
dir: folder of UNK files.
keep_bands: the band indexes to keep. Start from 1.
outdir: Defaults to 'truncated'.
"""
function truncate_unk(
    dir::String,
    keep_bands::AbstractVector{Int};
    outdir::String = "truncate",
)
    !isdir(outdir) && mkdir(outdir)

    regex = r"UNK(\d{5})\.\d"

    for unk in readdir(dir)
        match = match(regex, unk)
        match === nothing && continue

        println(unk)

        ik = parse(Int, match.captures[1])
        ik1, Ψ = read_unk(joinpath(dir, unk))
        @assert ik == ik1

        Ψ1 = Ψ[:, :, :, keep_bands]
        write_unk(joinpath(outdir, unk), ik, Ψ1)
    end

    nothing
end


"""
Truncate AMN/MMN/EIG/UNK(optional) files.

Args:
    seedname: seedname for input AMN/MMN/EIG files.
    keep_bands: Band indexes to be kept, start from 1.
    unk: Whether truncate UNK files. Defaults to false.
    outdir: output folder
"""
function truncate(
    seedname::String,
    keep_bands::AbstractVector{Int};
    unk::Bool = false,
    outdir::String = "truncate",
)
    @info "Truncat AMN/MMN/EIG files"

    !isdir(outdir) && mkdir(outdir)

    # E = read_eig("$seedname.eig")
    # n_bands = size(E, 1)
    # keep_bands = [i for i = 1:n_bands if i ∉ exclude_bands]

    truncate_mmn_eig(seedname, keep_bands, outdir)

    dir = dirname(seedname)

    if unk
        truncate_unk(dir, keep_bands, outdir)
    end

    println("Truncated files written in ", outdir)
    nothing
end
