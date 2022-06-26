import GarishPrint
import Configurations


Configurations.@option mutable struct InputParams
    # seedname.win
    seed_name::String = ""

    restart::String = ""

    read_amn::Bool = true

    # read the eig file (can be set to false when not disentangling)
    read_eig::Bool = true

    # write $file.optimize.amn at the end
    write_optimized_amn::Bool = true

    # will freeze if either n <= num_frozen or eigenvalue in frozen window
    dis_froz_num::Int = 0

    # 
    dis_froz_pao::Bool = false

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

    # Currently not read from toml, instead is used to save parameters from win
    kpath::Array{Float64,3} = zeros(3, 2, 1)
    kpath_label::Union{Matrix{Missing},Matrix{String}} = Matrix{Missing}(missing, 2, 1)
    bands_num_points::Int = 100
end

Base.show(io::IO, ::MIME"text/plain", x::InputParams) = GarishPrint.pprint_struct(io, x)

Base.@kwdef mutable struct Wannier90Nnkp
    # Note each column is a lattice vector, while in nnkp file each row is a lattice vector
    # 3 * 3
    real_lattice::Array{Float64,2}

    # 3 * 3
    recip_lattice::Array{Float64,2}

    n_kpts::Int

    # 3 * n_kpts
    kpoints::Array{Float64,2}

    # projections
    # auto_projections

    num_bvecs::Int

    # 4 * num_bvecs * n_kpts
    # nnkpts[1, ib, ik] = k+b equivalent kpoint k' in 1st BZ
    # nnkpts[2:4, ib, ik] = displacement vector from k' to k+b
    nnkpts::Array{Int,3}

    # exclude_bands
end

Base.@kwdef mutable struct InterpResults
    Obs_array::Array{ComplexF64,3}
    Obs_array_i::Array{ComplexF64,3}
    Obs_array_j::Array{ComplexF64,3}
    Uint::Array{ComplexF64,4}
    Uint_ik::Array{ComplexF64,4}
end


Base.@kwdef mutable struct Bands
    n_kpts::Int
    n_bands::Int

    # n_kpts, Cartesian, Å^-1 ?
    kpaths::Vector{Float64}
    # kpath coordinates, reduced, 3 * n_kpts
    kpaths_coord::Matrix{Float64}
    # n_kpts * n_bands
    energies::Matrix{Float64}

    # optional
    # index of high symmetry points
    num_symm_points::Int
    symm_points::Vector{Int}
    symm_points_label::Vector{String}

    # RGBA, 4 * n_kpts * n_bands
    # colors::Matrix{Int}
end

Base.@kwdef mutable struct AtomicWavefunction
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

Base.@kwdef mutable struct Projectabilities
    n_kpts::Int
    n_bands::Int
    # number of atomic wavefunctions
    num_wfcs::Int

    # atomic wavefunction types, size: num_wfcs
    wfcs_type::Vector{AtomicWavefunction}

    # projectability data, size: n_kpts * n_bands * num_wfcs
    proj::Array{Float64,3}
end
