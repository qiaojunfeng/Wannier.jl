using LinearAlgebra

export split_wannierize, split_unk, split_model

"""
    split_eig(E, U, n_val)

Split eigenvalues into two groups.

The separation is done by
1. construct Wannier gauge Hamiltonian,
    ``H_{\\bm{k}} = U_{\\bm{k}}^\\dagger [\\epsilon_{n \\bm{k}}] U_{\\bm{k}}``
2. diagonalize the Hamiltonian, the eigenvalues are sorted in ascending order,
    so that the first `n_val` eigenvalues are the occupied states,
    and the rest are the unoccupied states.

# Arguments
- `E`: eigenvalues
- `U`: (semi-)Unitary matrices gauge transformation

!!! tip

    Usually, the `U` are the maximally localized gauge matrices from a
    valence+conduction band Wannierization, this function split the
    gauge matrices `U` into two groups by diagonalizing the Hamiltonian,
    so we have two set of eigenvalues and gauge matrices for valence
    and conduction bands, respectively. However, the diagonalization
    introduce random gauges, so the returned two gauge matrices for
    valence and conduction are bad, we new to run a parallel transport
    to smoothen the gauges.
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
    split_unk(dir, Uv, Uc, outdir_val="val", outdir_cond="cond")

Rotate `UNK` files.

These are large matrices, so we read/write to disk for each kpoint sequentially,
inside the function.

# Arguments
- `dir`: directory where `UNK` files are stored
- `Uv`: the Wannier (semi-)unitary matrix for rotating valence bands
- `Uc`: the Wannier (semi-)unitary matrix for rotating conduction bands
- `outdir_val`: output directory for valence bands
- `outdir_cond`: output directory for conduction bands
"""
function split_unk(
    dir::AbstractString,
    Uv::Array{T,3},
    Uc::Array{T,3},
    outdir_val::AbstractString="val",
    outdir_cond::AbstractString="cond",
) where {T<:Complex}
    n_kpts = size(Uv, 3)
    n_kpts != size(Uc, 3) && error("incompatible n_kpts")
    n_val = size(Uv, 2)
    n_cond = size(Uc, 2)

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
        ΨUv = reshape(ΨUv, n_gx, n_gy, n_gz, n_val)
        ΨUc = reshape(ΨUc, n_gx, n_gy, n_gz, n_cond)

        val = joinpath(outdir_val, unk)
        write_unk(val, ik, ΨUv)

        cond = joinpath(outdir_cond, unk)
        write_unk(cond, ik, ΨUc)

        println("ik = ", ik, " files written: ", val, " ", cond)
    end

    return nothing
end

"""
    split_model(model, n_val)

Split the `Model` into two `Model`s.

# Arguments
- `model`: the `Model` to be split
- `n_val`: number of valence WFs

!!! note

    Rotation eigenvectors are also returned, useful for further
    rotation of `UNK` files or other operators.
"""
function split_model(model::Model, n_val::Int)
    n_wann = model.n_wann
    !(0 < n_val < n_wann) && error("n_val <= 0 or n_val >= n_wann")

    # EIG
    E = model.E
    U = model.A
    Ev, Ec, Vv, Vc = split_eig(E, U, n_val)
    UVv = rotate_A(U, Vv)
    UVc = rotate_A(U, Vc)

    model_v = rotate_gauge(model, UVv)
    @assert model_v.E ≈ Ev
    model_c = rotate_gauge(model, UVc)
    @assert model_c.E ≈ Ec

    return model_v, model_c, UVv, UVc
end

"""
    split_wannierize(model::Model, n_val::Int)

Split the model and run parallel transport to smoothen the gauge.

# Arguments
- `model`: the `Model` to be split
- `n_val`: number of valence WFs

!!! note

    Return two splitted `Model`s and rotation matrices which are
    useful for `UNK` files or other operators.
"""
function split_wannierize(model::Model, n_val::Int)
    model_v, model_c, UVv, UVc = split_model(model, n_val)

    Av, _ = parallel_transport(model_v)
    model_v.A .= Av
    Ac, _ = parallel_transport(model_c)
    model_c.A .= Ac

    return model_v, model_c, UVv, UVc
end
