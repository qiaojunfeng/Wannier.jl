using LinearAlgebra

export split_wannierize, split_unk, split_model

"""
    split_eig(E, U, eig_groups)

Split eigenvalues into several groups.

The separation is done by
1. construct Wannier gauge Hamiltonian,
    ``H_{\\bm{k}} = U_{\\bm{k}}^\\dagger [\\epsilon_{n \\bm{k}}] U_{\\bm{k}}``
2. diagonalize the Hamiltonian, the eigenvalues are sorted in ascending order,
    and they are split into several groups according to the indices in `eig_groups`.

# Arguments
- `E`: eigenvalues
- `U`: (semi-)Unitary matrices gauge transformation
- `eig_groups`: a Vector, in which each element is a Vector of indices of eigenvalues
    that belong to the same group.

!!! tip

    For example, if one wants to split valence+conduction into two groups,
    valence and conduction, respectively, then `eig_groups` should be
    `[[1:n_val], [n_val+1:n_wann]]`.
    Usually, the `U` are the maximally localized gauge matrices from a
    valence+conduction band Wannierization, this function split the
    gauge matrices `U` into two groups by diagonalizing the Hamiltonian,
    so we have two set of eigenvalues and gauge matrices for valence
    and conduction bands, respectively. However, the diagonalization
    introduce random gauges, so the returned two gauge matrices for
    valence and conduction are bad, we new to run a parallel transport
    to smoothen the gauges.
"""
function split_eig(
    E::Vector{Vector{T}}, U::Vector{Matrix{Complex{T}}}, eig_groups::AbstractVector{R}
) where {T<:Real,R<:AbstractVector{Int}}
    n_bands, n_wann = size(U[1])
    n_kpts = length(U)
    length(E[1]) != n_bands && error("incompatible n_bands")
    length(E) != n_kpts && error("incompatible n_kpts")

    n_groups = length(eig_groups)
    len_groups = [length(g) for g in eig_groups]
    n_wann == sum(len_groups) || error("incompatible eig_groups")
    E_groups = [[similar(E[1], len_groups[i]) for j in 1:n_kpts] for i in 1:n_groups]
    # Eigenvectors from diagonalization of Wannier Hamiltonian
    V_groups = [
        [similar(U[1], n_wann, len_groups[i]) for j in 1:n_kpts] for i in 1:n_groups
    ]

    # Since valence and conduction are split by a gap,
    # by diagonalizing the Wannier Hamiltonian and using the eigenvalues,
    # we can demix the WFs into two groups for valence and conduction, respectively.
    for ik in 1:n_kpts
        # Hamiltonian in WF basis
        Hₖ = Hermitian(U[ik]' * diagm(0 => E[ik]) * U[ik])

        # Diagonalize
        Dₖ, Vₖ = eigen(Hₖ)

        for ig in 1:n_groups
            E_groups[ig][ik] .= Dₖ[eig_groups[ig]]

            # Although the diagonalization destroy the smoothness of gauge,
            # one can run max localization to smooth these random gauge, since
            # there is no disentanglement. However, this requires many iterations
            # of max localization.
            # I store the gauge rotation due to diagonalization.
            V_groups[ig][ik] .= Vₖ[:, eig_groups[ig]]
        end
    end

    return collect(zip(E_groups, V_groups))
end

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
- `n_val`: number of valence states
"""
function split_eig(
    E::Vector{Vector{T}}, U::Vector{Matrix{Complex{T}}}, n_val::Int
) where {T<:Real}
    n_wann = size(U[1], 2)

    n_val < 1 && error("n_val < 0")
    n_val >= n_wann && error("n_val >= n_wann")

    eig_groups = [1:n_val, (n_val + 1):n_wann]
    return split_eig(E, U, eig_groups)
end

"""
    split_unk(dir, Us, outdirs)

Rotate `UNK` files.

These are large matrices, so we read/write to disk for each kpoint sequentially,
inside the function.

# Arguments
- `dir`: directory where `UNK` files are stored
- `Us`: a Vector of Wannier (semi-)unitary gauge matrices for rotating band groups
- `outdirs`: a Vector of output directories for each band group

# Keyword arguments
- `binary`: whether to write in Fortran binary format
"""
function split_unk(
    dir::AbstractString,
    Us::AbstractArray{T},
    outdirs::AbstractVector{R};
    binary::Bool=false,
) where {T<:AbstractArray{<:Complex,3},R<:AbstractString}
    length(Us) == length(outdirs) || error("incompatible Us and outdirs")
    n_kpts = size(Us[1], 3)
    all(length(U) == n_kpts for U in Us) || error("incompatible n_kpts")
    len_groups = [size(U[1], 2) for U in Us]

    println("UNK files will be written in: ")
    for (i, odir) in enumerate(outdirs)
        !isdir(odir) && mkdir(odir)
        @printf("    group %3d : %s", i, odir)
    end

    regex = r"UNK(\d{5})\.\d"

    for unk in readdir(dir)
        m = match(regex, unk)
        m === nothing && continue

        ik = parse(Int, m.captures[1])

        ik2, Ψ = read_unk(joinpath(dir, unk))
        @assert ik2 == ik

        # * does not support multiplication of high-dimensional arrays,
        # reshape it to matrix
        n_gx, n_gy, n_gz, n_bands = size(Ψ)
        Ψ = reshape(Ψ, :, n_bands)

        for (i, U) in enumerate(Us)
            Uₖ = U[ik]
            # rotate
            ΨU = Ψ * Uₖ
            # reshape back
            ΨU = reshape(ΨU, n_gx, n_gy, n_gz, len_groups[i])

            fname = joinpath(outdirs[i], unk)
            write_unk(fname, ik, ΨU; binary=binary)

            println("ik = ", ik, " files written: ", fname)
        end
    end

    return nothing
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

# Keyword arguments
- `binary`: whether to write in Fortran binary format
"""
function split_unk(
    dir::AbstractString,
    Uv::Vector{Matrix{T}},
    Uc::Vector{Matrix{T}},
    outdir_val::AbstractString="val",
    outdir_cond::AbstractString="cond";
    binary::Bool=false,
) where {T<:Complex}
    # !isdir(outdir_val) && mkdir(outdir_val)
    # !isdir(outdir_cond) && mkdir(outdir_cond)

    # println("UNK files will be written in: ")
    # println("    valence   : ", outdir_val)
    # println("    conduction: ", outdir_cond)

    return split_unk(dir, [Uv, Uc], [outdir_val, outdir_cond]; binary=binary)
end

"""
    split_model(model, eig_groups)

Split the `Model` into several `Model`s.

# Arguments
- `model`: the `Model` to be split
- `eig_groups`: a Vector, in which each element is a Vector of indices of eigenvalues
    that belong to the same group.

!!! note

    Rotation eigenvectors are also returned, useful for further
    rotation of `UNK` files or other operators.
"""
function split_model(
    model::Model, eig_groups::AbstractVector{R}
) where {R<:AbstractVector{Int}}
    E = model.E
    U = model.U
    EVs = split_eig(E, U, eig_groups)

    UVs = []
    models = []
    for EV in EVs
        UV = merge_gauge(U, EV[2])
        push!(UVs, UV)

        m = rotate_gauge(model, UV)
        @assert m.E ≈ EV[1]
        push!(models, m)
    end

    return collect(zip(models, UVs))
end

"""
    split_model(model, n_val)

Split the `Model` into two `Model`s.

# Arguments
- `model`: the `Model` to be split
- `n_val`: number of valence WFs
"""
function split_model(model::Model, n_val::Int)
    n_wann = model.n_wann
    !(0 < n_val < n_wann) && error("n_val <= 0 or n_val >= n_wann")

    eig_groups = [1:n_val, (n_val + 1):n_wann]
    return split_model(model, eig_groups)
end

"""
    split_wannierize(model::Model, eig_groups)

Split the model and run parallel transport to smoothen the gauge.

# Arguments
- `model`: the `Model` to be split
- `eig_groups`: a Vector, in which each element is a Vector of indices of eigenvalues
    that belong to the same group.

!!! note

    Return two separated `Model`s and rotation matrices which are
    useful for `UNK` files or other operators.
"""
function split_wannierize(
    model::Model, eig_groups::AbstractVector{R}
) where {R<:AbstractVector{Int}}
    model_Us = split_model(model, eig_groups)

    for (model, _) in model_Us
        U, _ = parallel_transport(model)
        model.U .= U
    end

    return model_Us
end

"""
    split_wannierize(model::Model, n_val::Int)

Split the model and run parallel transport to smoothen the gauge.

# Arguments
- `model`: the `Model` to be split
- `n_val`: number of valence WFs
"""
function split_wannierize(model::Model, n_val::Int)
    eig_groups = [1:n_val, (n_val + 1):(model.n_wann)]
    return split_wannierize(model, eig_groups)
end
