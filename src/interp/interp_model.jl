@doc raw"""
    struct InterpModel

The model for Wannier interpolation.

Store the real space matrices, e.g., the Hamiltonian ``H(\bm{R})``.

Usually after Wannierization of a [`Model`](@ref Model), we construct this
`InterpModel` for Wannier interpolation of operators.

# Fields
- `kRvectors`: the kpoints and R vectors
- `kpath`: the kpoint path for band structure
- `H`: `n_wann * n_wann * n_rvecs`, the Hamiltonian in real space
- `r`: `n_wann * n_wann * n_rvecs`, the position operator in real space
- `S`: `n_wann * n_wann * n_rvecs * 3`, optional, the spin operator

!!! note
    For MDRS interpolation, the `n_rvecs` in this docstring is actually the
    number of ``\tilde{\bm{R}}`` vectors, i.e., `KRVectors.Rvectors.n_r̃vecs`.
    For WS interpolation, it is just `KRVectors.Rvectors.n_rvecs`.
"""
struct InterpModel{T<:Real}
    # R vectors for Fourier transform
    kRvectors::KRVectors{T}

    # kpoint path for band structure
    kpath::KPath

    # Hamiltonian H(R), n_wann * n_wann * n_rvecs
    H::Array{Complex{T},3}

    # position operator r(R), n_wann * n_wann * n_rvecs * 3
    r::Array{Complex{T},4}

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

function InterpModel(
    kRvectors::KRVectors, kpath::KPath, H::Array{T,3}, r::Array{T,4}, S::Array{T,4}
) where {T<:Complex}
    n_wann = size(H, 1)
    return InterpModel(kRvectors, kpath, H, r, S, n_wann)
end

function _get_nrvecs(kRvectors::KRVectors)
    Rvecs = kRvectors.Rvectors
    if Rvecs isa RVectorsMDRS
        n_rvecs = Rvecs.n_r̃vecs
    elseif Rvecs isa RVectors
        n_rvecs = Rvecs.n_rvecs
    else
        error("Unknown Rvectors type: $(typeof(Rvecs))")
    end
    return n_rvecs
end

"""
    InterpModel(kRvectors, H, kpath, H, r)

A `InterpModel` constructor ignoring spin operator matrices.

# Arguments
- `kRvectors`: the kpoint and R vectors
- `kpath`: the kpoint path for band structure
- `H`: `n_wann * n_wann * n_rvecs`, the Hamiltonian in real space
- `r`: `n_wann * n_wann * n_rvecs * 3`, the position operator in real space
"""
function InterpModel(
    kRvectors::KRVectors, kpath::KPath, H::Array{T,3}, r::Array{T,4}
) where {T<:Complex}
    n_rvecs = _get_nrvecs(kRvectors)
    n_wann = size(H, 1)
    S = Array{T,4}(undef, n_wann, n_wann, n_rvecs, 3)

    return InterpModel(kRvectors, kpath, H, r, S)
end

"""
    InterpModel(kRvectors, H, kpath, H)

A `InterpModel` constructor ignoring position and spin operator matrices.

# Arguments
- `kRvectors`: the kpoint and R vectors
- `kpath`: the kpoint path for band structure
- `H`: `n_wann * n_wann * n_rvecs`, the Hamiltonian in real space
"""
function InterpModel(kRvectors::KRVectors, kpath::KPath, H::Array{T,3}) where {T<:Complex}
    n_rvecs = _get_nrvecs(kRvectors)
    n_wann = size(H, 1)
    r = Array{T,4}(undef, n_wann, n_wann, n_rvecs, 3)
    S = Array{T,4}(undef, n_wann, n_wann, n_rvecs, 3)
    return InterpModel(kRvectors, kpath, H, r, S)
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

    return InterpModel(kRvecs, kpath, Hᴿ)
end
