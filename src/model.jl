
# Base.@kwdef mutable struct Model
struct Model{T<:Real}
    # unit cell, 3 * 3, angstrom unit, each column is a lattice vector
    lattice::Mat3{T}

    # number of kpoints along 3 directions
    kgrid::Vec3{Int}

    # kpoints array, reduced coordinates, 3 * n_kpts
    # n_kpts is the last index since julia array is column-major
    kpoints::Matrix{T}

    # b vectors satisfying b1 condition
    bvectors::BVectors{T}

    # n_bands * n_kpts
    frozen_bands::BitMatrix

    # Mmn matrix, n_bands * n_bands * n_bvecs * n_kpts
    M::Array{Complex{T},4}

    # Amn matrix, n_bands * n_wann * n_kpts
    A::Array{Complex{T},3}

    # eigenvalues, n_bands * n_kpts
    E::Matrix{T}

    # spn matrix, n_bands * n_bands * 3 * n_kpts
    S::Array{Complex{T},4}

    # reciprocal cell, 3 * 3
    recip_lattice::Mat3{T}

    # number of bands
    n_bands::Int

    # number of Wannier functions (WFs)
    n_wann::Int

    # number of kpoints
    n_kpts::Int

    # number of neighbouring kpoints (i.e. b vectors)
    n_bvecs::Int
end

function Model(
    lattice::Mat3{T},
    kgrid::Vec3{Int},
    kpoints::Matrix{T},
    bvectors::BVectors{T},
    frozen_bands::AbstractMatrix{Bool},
    M::Array{Complex{T},4},
    A::Array{Complex{T},3},
    E::Matrix{T},
    S::Array{Complex{T},4},
) where {T<:Real}
    return Model(
        lattice,
        kgrid,
        kpoints,
        bvectors,
        BitMatrix(frozen_bands),
        M,
        A,
        E,
        S,
        get_recip_lattice(lattice),
        size(A, 1),
        size(A, 2),
        size(A, 3),
        bvectors.n_bvecs,
    )
end

function Model(
    lattice::Mat3{T},
    kgrid::Vec3{Int},
    kpoints::Matrix{T},
    bvectors::BVectors{T},
    frozen_bands::AbstractMatrix{Bool},
    M::Array{Complex{T},4},
    A::Array{Complex{T},3},
    E::Matrix{T},
) where {T<:Real}
    return Model(
        lattice,
        kgrid,
        kpoints,
        bvectors,
        frozen_bands,
        M,
        A,
        E,
        Array{Complex{T},4}(undef, 0, 0, 0, 0),
    )
end
