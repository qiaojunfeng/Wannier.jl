using Printf: @sprintf
using Dates: now

export write_chk

"""
    write_chk(filename, model, U; exclude_bands=nothing, binary=false)

Write the `model` to a Wannier90 `.chk` file, using the gauge `U`.

# Arguments
- `filename`: filename of the `.chk` file
- `model`: a `Model` struct
- `U`: `n_bands * n_wann * n_kpts` array for the gauge transformation

# Keyword Arguments
- `exclude_bands`: a list of band indices to exclude.
    This is irrelevant to the `model`, but `chk` file has this entry.
- `binary`: whether to write the `.chk` file in binary format.
"""
function write_chk(
    filename::AbstractString,
    model::Model,
    U::AbstractArray3;
    exclude_bands::Union{AbstractVector{Int},Nothing}=nothing,
    binary::Bool=false,
)
    header = @sprintf "Created by Wannier.jl %s" string(now())

    if isnothing(exclude_bands)
        exclude_bands = Vector{Int}()
    end

    checkpoint = "postwann"
    have_disentangled = true
    Ω = omega(model, U)
    dis_bands = trues(model.n_bands, model.n_kpts)

    # W90 has a special convention that the rotated Hamiltonian by the Uᵈ,
    # i.e., the Hamiltonian rotated by the gauge matrix from disentanglement,
    # needs to be diagonal. If I just store the whole unitary matrices Uᵈ*U
    # into the Uᵈ of chk, W90 will interpolate wrongly the band structure from
    # the written chk file.
    # Here, I first diagonalize the Hamiltonian Uᵈ'*E*Uᵈ, and store the
    # corresponding rotations into the U part of chk, to satisfy such convention.
    iszero(model.E) && error("E is all zero, cannot write chk file")

    H = get_Hk(model.E, U)
    if all(isdiag(H[:, :, ik]) for ik in axes(H, 3))
        # Uᵈ saved as disentanglement matrix, Uᵐ saved as max loc matrix
        Uᵈ = U
        Uᵐ = eyes_U(eltype(U), model.n_wann, model.n_kpts)
    else
        _, V = diag_Hk(H)
        Uᵈ = rotate_U(U, V)
        Uᵐ = permutedims(conj(V), [2, 1, 3])
    end

    M = rotate_M(model.M, model.bvectors.kpb_k, U)

    chk = WannierIO.Chk(
        header,
        exclude_bands,
        model.lattice,
        model.recip_lattice,
        model.kgrid,
        model.kpoints,
        checkpoint,
        have_disentangled,
        Ω.ΩI,
        dis_bands,
        Uᵈ,
        Uᵐ,
        M,
        Ω.r,
        Ω.ω,
    )

    return WannierIO.write_chk(filename, chk; binary)
end

"""
    write_chk(filename, model; exclude_bands=nothing, binary=false)

Write the `model` to a Wannier90 `.chk` file.

# Arguments
- `filename`: filename of the `.chk` file
- `model`: a `Model` struct

# Keyword Arguments
- `exclude_bands`: a list of band indices to exclude.
    This is irrelevant to the `model`, but `chk` file has this entry.
- `binary`: whether to write the `.chk` file in binary format.
"""
function write_chk(
    filename::AbstractString,
    model::Model;
    exclude_bands::Union{AbstractVector{Int},Nothing}=nothing,
    binary::Bool=false,
)
    return write_chk(filename, model, model.U; exclude_bands, binary)
end

"""
    Model(chk::Chk)

Construct a model from a `WannierIO.Chk` struct.

# Arguments
- `chk`: a `WannierIO.Chk` struct

# Keyword Arguments
- `E`: a `n_wann * n_kpts` array for the eigenvalues of the Hamiltonian.
    If not provided, it will be set to zero.
- `U`: a `n_wann * n_wann * n_kpts` array for the gauge transformation.

!!! warning

    The `Chk` struct does not contain eigenvalues, thus if `E` is not
    provided, it will be set to zero.

    Moreover, the `M` matrix in `Chk` is already rotated by the gauge
    transformation, thus by default, the `U` matrix is set to identity.
    Note that although maximal localization, or disentanglement (after
    frozen states are chosen), do not require eigenvalues (so the user
    can still Wannierize the `Model`), it is required when writing the
    `Model` to a `.chk` file, in [`write_chk`](@ref).

    Additionally, be careful that the `M` matrix is rotated, and this
    rotation needs to make sure that the rotated Hamiltonian is diagonal
    so that `E` stores the diagonal eigenvalues of the Hamiltonian.
"""
function Model(
    chk::WannierIO.Chk;
    E::Union{AbstractMatrix{Real},Nothing}=nothing,
    U::Union{AbstractArray3{Complex},Nothing}=nothing,
)
    atom_positions = zeros(Float64, 3, 0)
    atom_labels = Vector{String}()

    recip_lattice = get_recip_lattice(chk.lattice)
    # I try to generate bvectors, but it might happen that the generated bvectors
    # are different from the calculation corresponding to the chk file,
    # e.g. kmesh_tol is different
    bvectors = get_bvectors(chk.kpoints, recip_lattice)
    if bvectors.n_bvecs != chk.n_bvecs
        error("Number of bvectors is different from the number in the chk file")
    end

    frozen_bands = falses(chk.n_bands, chk.n_kpts)

    # the M in chk is already rotated by the U matrix
    M = chk.M
    # if nothing provided, I set U matrix as identity
    if isnothing(U)
        U = eyes_U(eltype(M), chk.n_wann, chk.n_kpts)
    end

    # no eig in chk file
    if isnothing(E)
        E = zeros(real(eltype(M)), chk.n_wann, chk.n_kpts)
    end

    model = Model(
        chk.lattice,
        atom_positions,
        atom_labels,
        chk.kgrid,
        chk.kpoints,
        bvectors,
        frozen_bands,
        M,
        U,
        E,
    )
    return model
end
