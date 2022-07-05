import GarishPrint
import Configurations


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
