export isisolated, isentangled

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
struct Model{T<:Real}
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
    Length-`n_kpoints` vector, each element is a length-`n_bvectors` vector, then
    each element is a `n_bands * n_bands` matrix.
    Also called `mmn` matrices in wannier90"""
    overlaps::Vector{Vector{Matrix{Complex{T}}}}

    """unitary or semi-unitary gauge transformation matrices, ``U_{\\bm{k}}``,
    or called the gauge matrices.
    Length-`n_kpoints` vector, each element is a `n_bands * n_wannier` matrix"""
    gauges::Vector{Matrix{Complex{T}}}

    """energy eigenvalues, ``\\varepsilon_{n \\bm{k}}``, length-`n_kpoints` vector,
    each element is a length-`n_bands` vector, in eV unit"""
    eigenvalues::Vector{Vector{T}}

    """mask for frozen bands. Length-`n_kpoints` vector, each element is a
    length-`n_bands` BitVector. If `true` the the state at that kpoint and band
    index is kept unchanged during the disentanglement procedure."""
    frozen_bands::Vector{BitVector}

    """mask for bands taking part in disentanglement. Length-`n_kpoints` vector,
    each element is a length-`n_bands` BitVector. If `true` the the state at that
    kpoint and band index participates the disentanglement procedure."""
    entangled_bands::Vector{BitVector}
end

# expose the fields of kstencil for convenience
function Base.propertynames(model::Model)
    return Tuple(
        [
            collect(fieldnames(typeof(model.kstencil)))
            collect(fieldnames(typeof(model)))
        ]
    )
end

function Base.getproperty(model::Model, sym::Symbol)
    type_stencil = typeof(getfield(model, :kstencil))
    if sym ∈ fieldnames(type_stencil)
        return getfield(getfield(model, :kstencil), sym)
    else
        # fallback
        return getfield(model, sym)
    end
end

function Model(
    lattice::AbstractMatrix,
    atom_positions::AbstractVector,
    atom_labels::AbstractVector,
    kstencil::KspaceStencil,
    overlaps::AbstractVector,
    gauges::AbstractVector,
    eigenvalues::AbstractVector,
    frozen_bands::AbstractVector,
    entangled_bands::AbstractVector=default_entangled_bands(gauges),
)
    # I introduce another set of accronyms for the dimensions of the matrices,
    # natoms -> n_atoms, nkpts -> n_kpoints, nbands -> n_bands, nwann -> n_wannier
    # so that they don't clash with those function names.
    natoms = length(atom_positions)
    @assert length(atom_labels) == natoms "atom_labels has wrong number of atoms"

    nkpts = n_kpoints(kstencil)
    @assert nkpts > 0 "empty kpoints"
    @assert length(overlaps) == nkpts "overlaps has wrong number of kpoints"
    @assert length(gauges) == nkpts "gauges has wrong number of kpoints"
    @assert length(eigenvalues) == nkpts "eigenvalues has wrong number of kpoints"
    @assert length(frozen_bands) == nkpts "frozen_bands has wrong number of kpoints"
    @assert length(entangled_bands) == nkpts "entangled_bands has wrong number of kpoints"
    @assert n_kpoints(kstencil) == nkpts "bvectors has wrong number of kpoints"

    nbvecs = n_bvectors(kstencil)
    @assert all(length(Mk) == nbvecs for Mk in overlaps) "overlaps has wrong number of b-vectors"

    nbands, nwann = size(gauges[1])
    @assert nbands ≥ nwann "number of bands must ≥ number of Wannier functions"
    @assert all(all(size(Mkb) == (nbands, nbands) for Mkb in Mk) for Mk in overlaps) "overlaps has wrong number of bands"
    @assert all(size(Uk) == (nbands, nwann) for Uk in gauges) "gauges has wrong number of bands or Wannier functions"
    @assert all(length(εk) == nbands for εk in eigenvalues) "eigenvalues has wrong number of bands"
    @assert all(length(fk) == nbands for fk in frozen_bands) "frozen_bands has wrong number of bands"
    @assert all(length(ek) == nbands for ek in entangled_bands) "entangled_bands has wrong number of bands"

    T = promote_type(eltype(lattice), eltype(recip_lattice))
    T = promote_type(T, eltype(real(overlaps[1][1])), eltype(real(gauges[1])))
    CT = Complex{T}

    return Model{T}(
        Mat3{T}(lattice),
        Vector{Vec3{T}}(atom_positions),
        atom_labels,
        kstencil,
        Vector{Vector{Matrix{CT}}}(overlaps),
        Vector{Matrix{CT}}(gauges),
        Vector{Vector{T}}(eigenvalues),
        frozen_bands,
        entangled_bands,
    )
end

n_atoms(model::Model) = length(model.atom_positions)
n_kpoints(model::Model) = n_kpoints(model.kstencil)
n_bvectors(model::Model) = n_bvectors(model.kstencil)
n_bands(model::Model) = isempty(model.gauges) ? 0 : size(model.gauges[1], 1)
n_wannier(model::Model) = isempty(model.gauges) ? 0 : size(model.gauges[1], 2)
real_lattice(model::Model) = model.lattice
reciprocal_lattice(model::Model) = reciprocal_lattice(model.kstencil)

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
    @printf(io, "  kgrid_size  =  %d %d %d\n", model.kgrid_size...)
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
