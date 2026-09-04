export BandEnergy

"""
    BandEnergy()

Request interpolated band energies in eV. The `band_energy` result field has
layout `n_wannier × n_kpoints`.
"""
struct BandEnergy <: Observable end

result_name(::BandEnergy) = :band_energy
requirements(::BandEnergy) = _ObservableRequirements(
    (:hamiltonian,), (:eigensystem,)
)

function _allocate_observable_result(
        ::BandEnergy,
        model::InterpolationModel,
        number_kpoints::Integer,
    )
    T = typeof(real(zero(eltype(model.operators.hamiltonian.coefficients))))
    return Matrix{T}(undef, n_wannier(model), number_kpoints)
end

function _assemble_observable!(
        destination::AbstractMatrix,
        ::BandEnergy,
        intermediates,
        workspace,
    )
    copyto!(destination, intermediates.eigenvalues)
    return destination
end
