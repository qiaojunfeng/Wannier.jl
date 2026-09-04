export BandEnergy, interpolate

"""Request interpolated band energies, in eV."""
struct BandEnergy end

function _band_energy!(destination::AbstractMatrix, hamiltonian::AbstractArray{<:Complex, 3})
    size(destination, 2) == size(hamiltonian, 3) ||
        throw(DimensionMismatch("destination and Hamiltonian batch sizes differ"))
    for batch_index in axes(hamiltonian, 3)
        matrix = Hermitian(view(hamiltonian, :, :, batch_index))
        destination[:, batch_index] .= eigvals(matrix)
    end
    return destination
end

"""
    interpolate(model, kpoints, BandEnergy())

Interpolate band energies at fractional-coordinate `kpoints`. The returned named
result has one field, `band_energy`, whose layout is
`n_wannier × n_kpoints` and whose unit is eV.
"""
function interpolate(
        model::InterpolationModel,
        kpoints::AbstractVector,
        ::BandEnergy,
    )
    number_kpoints = length(kpoints)
    number_kpoints > 0 || throw(ArgumentError("kpoints cannot be empty"))
    hamiltonian = model.operators.hamiltonian
    energy_type = typeof(real(zero(eltype(hamiltonian.coefficients))))
    band_energy = Matrix{energy_type}(undef, n_wannier(model), number_kpoints)

    batch_size = _interpolation_batch_size(model, number_kpoints)
    for first_index in 1:batch_size:number_kpoints
        kpoint_range = first_index:min(first_index + batch_size - 1, number_kpoints)
        hamiltonian_batch = _evaluate_real_space_operator(
            hamiltonian, model.real_space, kpoints, kpoint_range
        )
        _band_energy!(view(band_energy, :, kpoint_range), hamiltonian_batch)
    end
    return (; band_energy)
end

function interpolate(
        model::InterpolationModel,
        kpoints::AbstractVector,
        observables::Tuple,
    )
    length(observables) == 1 && only(observables) isa BandEnergy || throw(
        ArgumentError("this implementation slice supports only BandEnergy()"),
    )
    return interpolate(model, kpoints, only(observables))
end
