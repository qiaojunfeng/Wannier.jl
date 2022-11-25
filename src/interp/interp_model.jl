@doc raw"""
    struct InterpModel

The model for Wannier interpolation.

Store the real space matrices, e.g., the Hamiltonian ``H(\bm{R})``.

Usually after Wannierization of a [`Model`](@ref Model), we construct this
`InterpModel` for Wannier interpolation of operators.

# Fields
- `kRvectors`: the kpoints and R vectors
- `H`: `n_wann * n_wann * n_rvecs`, the Hamiltonian in real space
- `kpath`: the kpoint path for band structure
- `S`: `n_wann * n_wann * n_rvecs * 3`, optional, the spin operator
"""
struct InterpModel{T<:Real}
    # R vectors for Fourier transform
    kRvectors::KRVectors{T}

    # Hamiltonian H(R), n_wann * n_wann * n_rvecs
    H::Array{Complex{T},3}

    # kpoint path for band structure
    kpath::KPath

    # spin matrix S(R), n_wann * n_wann * n_rvecs * 3
    S::Array{Complex{T},4}

    # These fields are filled automatically upon construction

    # number of Wannier functions (WFs)
    n_wann::Int
end

function Base.getproperty(x::InterpModel, sym::Symbol)
    if sym ∈ fieldnames(KRVectors)
        return getfield(x.kRvectors, sym)
    elseif sym == :n_rvecs
        return getproperty(x.kRvectors, sym)
    else
        # fallback to getfield
        getfield(x, sym)
    end
end

function Base.show(io::IO, model::InterpModel)
    @info "n_wann = $(model.n_wann)"
    println(io, "")  # empty line

    @info "k & R vectors"
    show(io, model.kRvectors)
    println(io, "\n")

    @info "kpath"
    show(io, "text/plain", model.kpath)
    return nothing
end

function InterpModel(kRvectors, H, kpath, S)
    n_wann = size(H, 1)
    return InterpModel(kRvectors, H, kpath, S, n_wann)
end

"""
    InterpModel(kRvectors, H, kpath)

A `InterpModel` constructor ignoring spin operator matrices.

# Arguments
- `kRvectors`: the kpoint and R vectors
- `H`: `n_wann * n_wann * n_rvecs`, the Hamiltonian in real space
- `kpath`: the kpoint path for band structure
"""
function InterpModel(kRvectors, H, kpath)
    Rvecs = kRvectors.Rvectors
    if typeof(Rvecs) <: RVectorsMDRS
        n_rvecs = Rvecs.n_r̃vecs
    elseif typeof(Rvecs) <: RVectors
        n_rvecs = Rvecs.n_rvecs
    else
        error("Unknown Rvectors type: $(typeof(Rvecs))")
    end
    n_wann = size(H, 1)
    S = Array{eltype(H),4}(undef, n_wann, n_wann, n_rvecs, 3)
    return InterpModel(kRvectors, H, kpath, S)
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

    Hᵏ = get_Hk(model.E, model.A)
    Hᴿ = fourier(kRvecs, Hᵏ)

    return InterpModel(kRvecs, Hᴿ, kpath)
end
