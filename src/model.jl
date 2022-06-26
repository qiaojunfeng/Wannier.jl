
# Base.@kwdef mutable struct Model
struct Model
    # unit cell, 3 * 3, angstrom unit, each column is a lattice vector
    lattice::Mat3{Float64}

    # number of kpoints along 3 directions
    kgrid::Vec3{Int}

    # kpoints array, reduced coordinates, 3 * n_kpts
    # n_kpts is the last index since julia array is column-major
    kpoints::Matrix{Float64}

    # weights of each b vector, n_bvecs
    bvec_weights::Vector{Float64}

    # k+b vectors, k -> k + b (index of equivalent kpt in the 1st BZ), n_bvecs * n_kpts
    kpb_k::Matrix{Int}

    # displacements between k + b and k + b wrapped around into the recip_cell, 3 * n_bvecs * n_kpts
    kpb_b::Array{Int,3}

    # n_bands * n_kpts
    frozen_bands::BitMatrix

    # Mmn matrix, n_bands * n_bands * n_bvecs * n_kpts
    M::Array{ComplexF64,4}

    # Amn matrix, n_bands * n_wann * n_kpts
    A::Array{ComplexF64,3}

    # eigenvalues, n_bands * n_kpts
    E::Matrix{Float64}

    # spn matrix, n_bands * n_bands * 3 * n_kpts
    S::Array{ComplexF64,4}

    # reciprocal cell, 3 * 3
    recip_lattice::Mat3{Float64}

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
    lattice::Mat3{Float64},
    kgrid::Vec3{Int},
    kpoints::Matrix{Float64},
    bvec_weights::Vector{Float64},
    kpb_k::Matrix{Int},
    kpb_b::Array{Int,3},
    frozen_bands::BitMatrix,
    M::Array{ComplexF64,4},
    A::Array{ComplexF64,3},
    E::Matrix{Float64},
    S::Array{ComplexF64,4},
)
    Model(
        lattice = lattice,
        kgrid = kgrid,
        kpoints = kpoints,
        bvec_weights = bvec_weights,
        kpb_k = kpb_k,
        kpb_b = kpb_b,
        frozen_bands = frozen_bands,
        M = M,
        A = A,
        E = E,
        S = S,
        recip_lattice = gett_recip_lattice(lattice),
        n_bands = size(A, 1),
        n_wann = size(A, 2),
        n_kpts = size(A, 3),
        n_bvecs = length(bvec_weights),
    )
end
