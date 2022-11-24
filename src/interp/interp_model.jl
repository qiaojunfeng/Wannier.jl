"""
    struct InterpModel

The model for Wannier interpolation.

# Fields
- `model`: the Wannierization [`Model`](@ref Model)
- `kRvectors`: the kpoint and R vectors
- `kpath`: the kpoint path for band structure
- `S`: `n_bands * n_bands * n_kpts * 3`, the spin operator matrices
"""
struct InterpModel{T<:Real}
    # model storing results of Wannierization
    model::Model{T}

    # R vectors for Fourier transform
    kRvectors::KRVectors{T}

    # kpoint path for band structure
    kpath::KPath

    # spin matrix, n_bands * n_bands * n_kpts * 3
    S::Array{Complex{T},4}
end

function Base.getproperty(x::InterpModel, sym::Symbol)
    if sym ∈ fieldnames(Model)
        return getfield(x.model, sym)
    elseif sym ∈ fieldnames(KRVectors)
        return getfield(x.kRvectors, sym)
    elseif sym == :n_rvecs
        return getproperty(x.kRvectors, sym)
    else
        # fallback to getfield
        getfield(x, sym)
    end
end

function Base.show(io::IO, model::InterpModel)
    @info "model"
    show(io, model.model)
    println(io, "\n")

    @info "k & R vectors"
    show(io, model.kRvectors)
    println(io, "\n")

    @info "kpath"
    show(io, "text/plain", model.kpath)
    return nothing
end

"""
    InterpModel(
        model::Model{T}, kRvectors::KRVectors{T}, kpath::KPath
    ) where {T<:Real}

A `InterpModel` constructor ignoring spin operator matrices.

# Arguments
- `model`: the Wannierization [`Model`](@ref Model)
- `kRvectors`: the kpoint and R vectors
- `kpath`: the kpoint path for band structure
"""
function InterpModel(model::Model{T}, kRvectors::KRVectors{T}, kpath::KPath) where {T<:Real}
    n_bands = model.n_bands
    n_kpts = model.n_kpts
    return InterpModel(
        model, kRvectors, kpath, Array{Complex{T},4}(undef, n_bands, n_bands, 3, n_kpts)
    )
end

"""
    InterpModel(model::Model; mdrs::Bool=true)

Construct a `InterpModel` from a Wannierization [`Model`](@ref Model).

The `kpath` will be auto generated from the lattice by using
[`get_kpath`](@ref get_kpath).

# Arguments
- `model`: the Wannierization [`Model`](@ref Model)

# Keyword Arguments
- `mdrs`: whether to use MDRS interpolation
"""
function InterpModel(model::Model; mdrs::Bool=true)
    if mdrs
        centers = center(model)
        # from cartesian to fractional
        centers = inv(model.lattice) * centers
        Rvecs = get_Rvectors_mdrs(model.lattice, model.kgrid, centers)
    else
        Rvecs = get_Rvectors_ws(model.lattice, model.kgrid)
    end
    kRvecs = KRVectors(model.lattice, model.kgrid, model.kpoints, Rvecs)

    kpath = get_kpath(model.lattice, model.atom_positions, model.atom_labels)

    return InterpModel(model, kRvecs, kpath)
end
