# to extend these methods
import WannierIO: write_chk

"""
    $(SIGNATURES)

Write a `Model` to a wannier90 `chk` file.

# Arguments
- `filename`: filename of the `.chk` file
- `model`: a `Model` struct
- `U`: gauge transformation matrices, default is `model.U`.

# Keyword Arguments
- `exclude_bands`: a list of band indices to be excluded.
- `binary`: write the `chk` file in binary format
- `header`: header of the `chk` file

!!! note

    The `exclude_bands` is totally not used throughout the code.
    However, if one want to use wannier90 to restart from the written `chk` file,
    the `exclude_bands` must be consistent with that in wannier90 `win` file.
    The wannier90 `exclude_bands` input parameter in only used in the pw2wannier90
    step to remove some bands (usually semicore states) when computing `amn`/`mmn`/`eig`
    files, therefore this `exclude_bands` is totally irrelevant to our `Model` struct.
"""
function write_chk(
    filename::AbstractString,
    model::Model,
    gauges::AbstractVector=model.gauges;
    exclude_bands::AbstractVector=default_exclude_bands(),
    binary::Bool=false,
    header=default_header(),
)
    checkpoint = "postwann"
    have_disentangled = isentangled(model)
    Ω = omega(model, gauges)

    if have_disentangled
        # W90 has a special convention that the rotated Hamiltonian (by the Udis),
        # i.e., the Hamiltonian rotated by the gauge matrix from disentanglement,
        # needs to be diagonal. If I just store the whole unitary matrices Uᵈ*Uml
        # into the Udis of chk, W90 will interpolate wrongly the band structure from
        # the written chk file.
        # Here, I first diagonalize the Hamiltonian Udis' * E * Udis, and store the
        # corresponding rotations into the Uml part of chk, to satisfy such convention.
        @assert !iszero(model.eigenvalues) "E is all zero, cannot write chk file"

        H = transform_gauge(model.eigenvalues, gauges)
        if all(h -> isdiag(h), H)
            # Udis saved as disentanglement matrix, Uml saved as max loc matrix
            Udis = gauges
            Uml = identity_gauge(eltype(gauges), n_kpoints(model), n_wannier(model))
        else
            V = map(h -> eigen(h).vectors, H)
            Udis = merge_gauge(gauges, V)
            Uml = map(v -> v', V)
        end
    else
        Udis = zero_gauge(eltype(gauges), n_kpoints(model), n_wannier(model))
        Uml = gauges
    end
    overlaps = transform_gauge(model.overlaps, model.kstencil.kpb_k, gauges)

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
        model.dis_bands,
        Udis,
        Uml,
        overlaps,
        Ω.r,
        Ω.ω,
    )

    return WannierIO.write_chk(filename, chk; binary)
end

"""
    $(SIGNATURES)

Construct a `Model` from a `WannierIO.Chk` struct.

# Arguments
- `chk`: a `WannierIO.Chk` struct

!!! warning

    The `Chk` struct does not contain eigenvalues, and the `Model.E`
    will be set to zeros.

    Moreover, the `M` matrix in `Chk` is already rotated by the gauge
    transformation, thus by default, the `U` matrix is set to identity.
    Note that although maximal localization, or disentanglement (after
    frozen states are chosen), do not require eigenvalues (so the user
    can still Wannierize the `Model`), it is required when writing the
    `Model` to a `chk` file, in [`write_chk`](@ref).

    Additionally, be careful that the `M` matrix is rotated, and this
    rotation needs to make sure that the rotated Hamiltonian is diagonal
    so that `E` stores the diagonal eigenvalues of the Hamiltonian.
"""
function Model(chk::WannierIO.Chk; kmesh_tol=default_w90_kmesh_tol())
    atom_positions = Vec3{Float64}[]
    atom_labels = Vector{String}()
    @warn "chk file does not contain info on atom positions and labels, set them to empty"

    # I try to generate bvectors, but it might happen that the generated bvectors
    # are different from the calculation corresponding to the chk file,
    # e.g. kmesh_tol is different
    kgrid = KpointGrid(chk.recip_lattice, chk.kgrid, chk.kpoints)
    kstencil = generate_stencil(kgrid; atol=kmesh_tol)
    @warn "The generated bvectors might be different from that used in the " *
        "chk file, if the wannier90 input parameter `kmesh_tol` is different from " *
        "its default value."
    @assert n_bvectors(kstencil) == chk.n_bvecs "n_bvecs different from chk file"

    frozen_bands = [falses(chk.n_bands) for _ in 1:(chk.n_kpts)]
    @warn "chk file does not contain info on frozen bands, set all to false"

    # the M in chk is already rotated by the U matrix
    overlaps = chk.M
    # set U matrix as identity
    T = eltype(overlaps[1][1])
    gauges = identity_gauge(T, chk.n_kpts, chk.n_wann)
    @warn "chk file only contains `overlaps` matrix rotated by the gauge " *
        "transformation, thus the `gauges` matrix is set to identity."

    # no eig in chk file
    eigenvalues = [zeros(real(T), chk.n_wann) for _ in 1:(chk.n_kpts)]
    @warn "chk file does not contain energy eigenvalues, set all to zero"

    return Model(
        chk.lattice,
        atom_positions,
        atom_labels,
        kgrid,
        kstencil,
        overlaps,
        gauges,
        eigenvalues,
        frozen_bands,
    )
end
