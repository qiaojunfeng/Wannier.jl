export isisolated, isentangled
export kpoints, kgrid_size, kpb_k, kpb_G, bweights

"""
    $(TYPEDEF)

A high-level data structure containing necessary parameters and matrices for
constructing Wannier functions (WFs), or called, Wannierization.

In general, the problem of Wannierization is to find a set of unitary matrices
``U_{\\bm{k}}`` that gives a localized representation of the Bloch states
``\\psi_{n \\bm{k}}``. Depending on the inputs, the Wannierization problem can be
categorized into two classes:
- isolated manifold: when number of bands = number of Wannier functions
- entangled manifold: when number of bands > number of Wannier functions

# Fields
Using these accronyms,
    - `n_atoms`: number of atoms
    - `n_kpoints`: number of kpoints
    - `n_bvectors`: number of b-vectors
    - `n_bands`: number of bands
    - `n_wannier`: number of Wannier functions
the fields are defined as follows:
$(FIELDS)
"""
struct Model{T <: Real}
    """unit cell, 3 * 3, each column is a lattice vector in Å unit"""
    lattice::Mat3{T}

    """atomic positions, length-`n_atoms` vector, each element is a `Vec3` of
    fractional coordinates"""
    atom_positions::Vector{Vec3{T}}

    """atomic labels, length-`n_atoms` vector of string"""
    atom_labels::Vector{String}

    """stencil for finite differences on the kpoint grid, also called
    ``\\mathbf{b}``-vectors. Should satisfy completeness condition, see
    [`KspaceStencil`](@ref)"""
    kstencil::KspaceStencil{T}

    """overlap matrices between neighboring wavefunctions, ``M_{\\bm{k},\\bm{b}}``.
    `n_bands × n_bands × n_bvectors × n_kpoints` array.
    Also called `mmn` matrices in wannier90"""
    overlaps::Array{Complex{T}, 4}

    """unitary or semi-unitary gauge transformation matrices, ``U_{\\bm{k}}``,
    or called the gauge matrices.
    `n_bands × n_wannier × n_kpoints` array"""
    gauges::Array{Complex{T}, 3}

    """energy eigenvalues, ``\\varepsilon_{n \\bm{k}}``,
    `n_bands × n_kpoints` matrix, in eV unit"""
    eigenvalues::Matrix{T}

    """mask for frozen bands. `n_bands × n_kpoints` BitMatrix.
    If `true` the the state at that kpoint and band index is kept unchanged
    during the disentanglement procedure."""
    frozen_bands::BitMatrix

    """mask for bands taking part in disentanglement. `n_bands × n_kpoints` BitMatrix.
    If `true` the the state at that kpoint and band index participates the
    disentanglement procedure."""
    entangled_bands::BitMatrix
end

function Model(
        lattice::AbstractMatrix,
        atom_positions::AbstractVector,
        atom_labels::AbstractVector,
        kstencil::KspaceStencil,
        overlaps::AbstractArray{<:Complex, 4},
        gauges::AbstractArray{<:Complex, 3},
        eigenvalues::AbstractMatrix{<:Real},
        frozen_bands::AbstractMatrix{Bool},
        entangled_bands::AbstractMatrix{Bool} = default_entangled_bands(gauges),
    )
    natoms = length(atom_positions)
    @assert length(atom_labels) == natoms "atom_labels has wrong number of atoms"

    nkpts = n_kpoints(kstencil)
    nbvecs = n_bvectors(kstencil)
    @assert nkpts > 0 "empty kpoints"
    @assert size(overlaps, 4) == nkpts "overlaps has wrong number of kpoints"
    @assert size(overlaps, 3) == nbvecs "overlaps has wrong number of b-vectors"
    @assert size(gauges, 3) == nkpts "gauges has wrong number of kpoints"
    @assert size(eigenvalues, 2) == nkpts "eigenvalues has wrong number of kpoints"
    @assert size(frozen_bands, 2) == nkpts "frozen_bands has wrong number of kpoints"
    @assert size(entangled_bands, 2) == nkpts "entangled_bands has wrong number of kpoints"

    nbands = size(overlaps, 1)
    nwann = size(gauges, 2)
    @assert size(overlaps, 2) == nbands "overlaps must be square in band dimensions"
    @assert nbands ≥ nwann "number of bands must ≥ number of Wannier functions"
    @assert size(gauges, 1) == nbands "gauges has wrong number of bands"
    @assert size(eigenvalues, 1) == nbands "eigenvalues has wrong number of bands"
    @assert size(frozen_bands, 1) == nbands "frozen_bands has wrong number of bands"
    @assert size(entangled_bands, 1) == nbands "entangled_bands has wrong number of bands"

    T = promote_type(eltype(lattice), real(eltype(overlaps)), real(eltype(gauges)))
    CT = Complex{T}

    return Model{T}(
        Mat3{T}(lattice),
        Vector{Vec3{T}}(atom_positions),
        atom_labels,
        kstencil,
        Array{CT, 4}(overlaps),
        Array{CT, 3}(gauges),
        Matrix{T}(eigenvalues),
        BitMatrix(frozen_bands),
        BitMatrix(entangled_bands),
    )
end

"""
    $(SIGNATURES)

Construct a `Model` from an existing `Model`, reuse the lattice and atomic
information, but change the `kstencil` and `overlaps`.

For instance, use a cubic-6-neighbors `KspaceStencil` for
[`parallel_transport`](@ref).
"""
function Model(
        model::Model,
        kstencil::KspaceStencil,
        overlaps::AbstractArray{<:Complex, 4},
        gauges::AbstractArray{<:Complex, 3} = model.gauges,
        eigenvalues::AbstractMatrix{<:Real} = model.eigenvalues,
        frozen_bands::AbstractMatrix{Bool} = model.frozen_bands,
    )
    return Model(
        model.lattice,
        model.atom_positions,
        model.atom_labels,
        kstencil,
        overlaps,
        gauges,
        eigenvalues,
        frozen_bands,
    )
end

n_atoms(model::Model) = length(model.atom_positions)
n_kpoints(model::Model) = n_kpoints(model.kstencil)
n_bvectors(model::Model) = n_bvectors(model.kstencil)
n_bands(model::Model) = size(model.gauges, 1)
n_wannier(model::Model) = size(model.gauges, 2)
CrystalBase.real_lattice(model::Model) = model.lattice
CrystalBase.reciprocal_lattice(model::Model) = reciprocal_lattice(model.kstencil)

"""
    kpoints(::Model)
    kpoints(::KspaceStencil)

Fractional coordinates of the kpoint grid.
"""
kpoints(model::Model) = model.kstencil.kpoints
kpoints(kstencil::KspaceStencil) = kstencil.kpoints

"""
    kgrid_size(::Model)
    kgrid_size(::KspaceStencil)

Number of kpoints along the three reciprocal lattice vectors.
"""
kgrid_size(model::Model) = model.kstencil.kgrid_size
kgrid_size(kstencil::KspaceStencil) = kstencil.kgrid_size

"""
    kpb_k(::Model)
    kpb_k(::KspaceStencil)

Indices of ``\\mathbf{k+b}`` kpoints.
"""
kpb_k(model::Model) = model.kstencil.kpb_k
kpb_k(kstencil::KspaceStencil) = kstencil.kpb_k

"""
    kpb_G(::Model)
    kpb_G(::KspaceStencil)

Reciprocal lattice displacements for ``\\mathbf{k+b}`` kpoints.
"""
kpb_G(model::Model) = model.kstencil.kpb_G
kpb_G(kstencil::KspaceStencil) = kstencil.kpb_G

"""
    bweights(::Model)
    bweights(::KspaceStencil)

Weights of the ``\\mathbf{b}``-vectors.
"""
bweights(model::Model) = model.kstencil.bweights
bweights(kstencil::KspaceStencil) = kstencil.bweights

"""
    $(SIGNATURES)

Is entangled manifold?
"""
isentangled(model::Model) = n_bands(model) > n_wannier(model)

"""
    $(SIGNATURES)

Is isolated manifold?
"""
isisolated(model::Model) = n_bands(model) == n_wannier(model)

function Base.show(io::IO, ::MIME"text/plain", model::Model)
    show_lattice(io, model.lattice)
    println(io)

    @printf(io, "atoms:             fractional\n")
    for (i, (label, pos)) in enumerate(zip(model.atom_labels, model.atom_positions))
        @printf(io, " %3d  %3s: %9.6f %9.6f %9.6f\n", i, label, pos...)
    end
    println(io, repeat("-", 80))

    show(io, MIME"text/plain"(), model.kstencil)
    println(io)
    println(io, repeat("-", 80))

    println(io, "Summary:")
    @printf(io, "  kgrid_size  =  %d %d %d\n", kgrid_size(model)...)
    @printf(io, "  n_kpoints   =  %d\n", n_kpoints(model))
    @printf(io, "  n_bvectors  =  %d\n", n_bvectors(model))
    @printf(io, "  n_bands     =  %d\n", n_bands(model))
    @printf(io, "  n_wannier   =  %d", n_wannier(model))
    return nothing
end

"""
    $(SIGNATURES)

Compare two `Model` objects.
"""
function Base.isapprox(a::Model, b::Model; kwargs...)
    return isapprox_struct(a, b; kwargs...)
end
