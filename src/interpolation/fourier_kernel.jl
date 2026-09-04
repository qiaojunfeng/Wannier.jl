const _DEFAULT_INTERPOLATION_MEMORY = 64 * 1024 * 1024
const _MAXIMUM_INTERPOLATION_BATCH = 4096

function _interpolation_batch_size(model::InterpolationModel, number_kpoints::Integer)
    number_kpoints > 0 || return 0
    T = eltype(model.operators.hamiltonian.coefficients)
    bytes_per_kpoint = sizeof(T) * (n_Rvectors(model) + n_wannier(model)^2)
    memory_limited = max(1, _DEFAULT_INTERPOLATION_MEMORY ÷ max(1, bytes_per_kpoint))
    return min(number_kpoints, memory_limited, _MAXIMUM_INTERPOLATION_BATCH)
end

function _fourier_phase_block(
        real_space::RealSpaceDomain,
        kpoints::AbstractVector,
        kpoint_range::AbstractUnitRange,
        ::Type{T},
    ) where {T <: Complex}
    phase = Matrix{T}(undef, length(real_space), length(kpoint_range))
    for (batch_index, kpoint_index) in enumerate(kpoint_range)
        kpoint = kpoints[kpoint_index]
        length(kpoint) == 3 || throw(ArgumentError("each k point must have length three"))
        for (vector_index, vector) in enumerate(real_space.vectors)
            phase[vector_index, batch_index] = cis(2π * dot(kpoint, vector))
        end
    end
    return phase
end

function _evaluate_real_space_operator(
        operator::RealSpaceOperator,
        real_space::RealSpaceDomain,
        kpoints::AbstractVector,
        kpoint_range::AbstractUnitRange,
    )
    coefficients = operator.coefficients
    number_vectors = size(coefficients, ndims(coefficients))
    number_vectors == length(real_space) ||
        throw(DimensionMismatch("operator and real-space domains differ"))

    T = eltype(coefficients)
    phase = _fourier_phase_block(real_space, kpoints, kpoint_range, T)
    packed_coefficients = reshape(coefficients, :, number_vectors)
    packed_values = Matrix{T}(
        undef, size(packed_coefficients, 1), length(kpoint_range)
    )
    mul!(packed_values, packed_coefficients, phase)
    return reshape(
        packed_values, size(coefficients)[1:(end - 1)]..., length(kpoint_range)
    )
end
