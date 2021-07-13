module Parameters

import Configurations

mutable struct CoreData
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

    # num_bands * num_kpts
    frozen::BitMatrix

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

Configurations.@option struct InputParams
    # seedname.win
    seed_name::String = ""

    read_amn::Bool = true

    # read the eig file (can be set to false when not disentangling)
    read_eig::Bool = true

    # write $file.optimize.amn at the end
    write_optimized_amn::Bool = true

    # will freeze if either n <= num_frozen or eigenvalue in frozen window
    dis_froz_num::Int = 0

    # dis_froz_min <= energy <= dis_froz_max bands are frozen
    dis_froz_win::Bool = false
    dis_froz_win_min::Float64 = -Inf
    dis_froz_win_max::Float64 = -Inf

    # select frozen bands based on projectability
    dis_froz_proj::Bool = false
    # TODO: < threshold bands are excluded
    dis_froz_proj_min::Float64 = 0.1
    # >= threshold bands are frozen
    dis_froz_proj_max::Float64 = 0.9

    # will also freeze additional eigenvalues if the freezing cuts a cluster. Set to 0 to disable
    dis_froz_degen::Bool = false
    dis_froz_degen_thres::Float64 = 1e-6

    fix_centers::Bool = false
    # TOML only accepts Vector of Vector for matrix input,
    # so each row is a coord vector, however almost all the internal matrices are column vector
    fix_centers_coords::Vector{Vector{Float64}} = []
    # [ # angstrom
    # 1.34940   1.34940   1.34940; 
    # 1.34940   1.34940   1.34940; 
    # 1.34940   1.34940   1.34940;
    # 1.34940   1.34940   1.34940;
    # 0.00000   0.00000   0.00000;
    # 0.00000   0.00000   0.00000;
    # 0.00000   0.00000   0.00000;
    # 0.00000   0.00000   0.00000;
    # ]

    # expert/experimental features
    # perform a global rotation by a phase factor at the end
    do_normalize_phase::Bool = false

    # randomize initial gauge
    do_randomize_gauge::Bool = false

    # only minimize sum_n <r^2>_n, not sum_n <r^2>_n - <r>_n^2
    only_r2::Bool = false

    # tolerance on spread
    omega_tol::Float64 = 1e-10
    # tolerance on gradient
    gradient_tol::Float64 = 1e-4
    # maximum optimization iterations
    max_iter::Int = 150 # 3000
    # history size of BFGS
    history_size::Int = 100
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
