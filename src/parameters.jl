
module Parameters

using Base: String
mutable struct InputParams
    # seedname.win
    seed_name::String

    # unit cell, 3 * 3, angstrom unit, each column is a lattice vector
    unit_cell::Array{Float64,2}
    # reciprocal cell, 3 * 3
    recip_cell::Array{Float64,2}

    # number of bands
    num_bands::Int
    # number of Wannier functions (WFs)
    num_wann::Int

    # number of kpoints
    num_kpts::Int
    # number of kpoints along 3 directions
    kpts_size::Array{Int,1}
    # kpoints array, reduced coordinates, 3 * num_kpts
    # num_kpts is the last index since julia array is column-major
    kpts::Array{Float64,2}

    # number of neighbouring kpoints (i.e. b vectors)
    num_bvecs::Int
    # k+b vectors, k -> k + b (index of equivalent kpt in the 1st BZ), nbvecs * num_kpts
    kpbs::Array{Int,2}
    # weights of each b vector, nbvecs * num_kpts
    kpbs_weight::Array{Float64,2}
    # displacements between k + b and k + b wrapped around into the recipcell, 3 * nbvecs * num_kpts
    kpbs_disp::Array{Int,3}

    # Mmn matrix, nbands * nbands * nbvecs * num_kpts
    mmn::Array{ComplexF64,4}
    # Amn matrix, nbands * nwann * num_kpts
    amn::Array{ComplexF64,3}
    # eigenvalues, nbands * num_kpts
    eig::Array{Float64,2}
    # spn matrix, nbands * nbands * 3 * num_kpts
    # spn::Array{ComplexF64,4}

    # TODO: remove this?
    # map::Bool
    # logMethod::Bool
end

mutable struct Wannier90Nnkp
    # Note each column is a lattice vector, while in nnkp file each row is a lattice vector
    # 3 * 3
    real_lattice::Array{Float64,2}

    # 3 * 3
    recip_lattice::Array{Float64,2}

    num_kpts::Int

    # 3 * num_kpts
    kpoints::Array{Float64,2}

    # projections
    # auto_projections

    num_bvecs::Int

    # 4 * num_bvecs * num_kpts
    # nnkpts[1, ib, ik] = k+b equivalent kpoint k' in 1st BZ
    # nnkpts[2:4, ib, ik] = displacement vector from k' to k+b
    nnkpts::Array{Int,3}

    # exclude_bands
end

mutable struct InterpResults
    Obs_array::Array{ComplexF64,3}
    Obs_array_i::Array{ComplexF64,3}
    Obs_array_j::Array{ComplexF64,3}
    Uint::Array{ComplexF64,4}
    Uint_ik::Array{ComplexF64,4}
end


mutable struct Bands
    num_kpts::Int
    num_bands::Int

    # num_kpts
    kpaths::Vector{Float64}
    # kpath coordinates, reduced, 3 * num_kpts
    kpaths_coord::Matrix{Float64}
    # num_kpts * num_bands
    energies::Matrix{Float64}

    # optional
    # index of high symmetry points
    num_symm_points::Int
    symm_points::Vector{Int}
    symm_points_label::Vector{String}

    # RGBA, 4 * num_kpts * num_bands
    # colors::Matrix{Int}
end

mutable struct AtomicWavefunction
    # to differentiate same kind of atoms
    atom_index::Int
    # e.g. "Si", use "Si1", "Si2" to differentiate same kind but different types (spin up, down)
    atom_label::String
    # orbital label, e.g. "3S"
    wfc_label::String
    
    # quantum numbers
    n::Int
    l::Int
    m::Int
end

mutable struct Projectabilities
    num_kpts::Int
    num_bands::Int
    # number of atomic wavefunctions
    num_wfcs::Int

    # atomic wavefunction types, size: num_wfcs
    wfcs_type::Vector{AtomicWavefunction}

    # projectability data, size: num_kpts * num_bands * num_wfcs
    proj::Array{Float64,3}
end

end
