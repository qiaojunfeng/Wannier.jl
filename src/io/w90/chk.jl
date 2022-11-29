using Printf: @sprintf
using Dates: now

export write_chk

"""
    write_chk(filename, model; exclude_bands=nothing)

Write the `model` to a Wannier90 `.chk` file.

# Arguments
- `filename`: filename of the `.chk` file
- `model`: a `Model` struct

# Keyword Arguments
- `exclude_bands`: a list of band indices to exclude.
    This is irrelevant to the `model`, but `chk` file has this entry.
"""
function write_chk(
    filename::AbstractString,
    model::Model;
    exclude_bands::Union{AbstractVector{Int},Nothing}=nothing,
    binary::Bool=false,
)
    header = @sprintf "Created by Wannier.jl %s" string(now())

    if isnothing(exclude_bands)
        exclude_bands = Vector{Int}()
    end

    checkpoint = "postwann"
    have_disentangled = true
    Ω = omega(model)
    dis_bands = trues(model.n_bands, model.n_kpts)
    Uᵈ = model.U
    U = eyes_U(eltype(Uᵈ), model.n_wann, model.n_kpts)
    M = rotate_M(model.M, model.bvectors.kpb_k, model.U)

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
        U,
        M,
        Ω.r,
        Ω.ω,
    )

    WannierIO.write_chk(filename, chk; binary=binary)
    return nothing
end

"""
    Model(chk::Chk)

Construct a model from a `WannierIO.Chk` struct.
"""
function Model(chk::WannierIO.Chk)
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
    # so I set U matrix as identity
    U = eyes_U(eltype(M), chk.n_wann, chk.n_kpts)

    # no eig in chk file
    E = zeros(Float64, chk.n_wann, chk.n_kpts)

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
