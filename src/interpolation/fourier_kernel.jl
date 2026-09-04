const _DEFAULT_INTERPOLATION_MEMORY = 64 * 1024 * 1024
const _MAXIMUM_INTERPOLATION_BATCH = 4096

function _interpolation_batch_size(
        model::InterpolationModel,
        number_kpoints::Integer;
        operator_names::Tuple = (:hamiltonian,),
        derivative_operator_names::Tuple = (),
    )
    number_kpoints > 0 || return 0
    T = eltype(model.operators.hamiltonian.coefficients)
    number_wannier = n_wannier(model)
    number_vectors = n_Rvectors(model)
    number_complex_values = number_vectors + 2number_wannier^2
    if !isempty(derivative_operator_names)
        number_complex_values += number_vectors
        for operator_name in derivative_operator_names
            operator = getproperty(model.operators, operator_name)
            number_complex_values +=
                3prod(size(operator.coefficients)[1:(end - 1)])
        end
    end
    for operator_name in operator_names
        operator_name == :hamiltonian && continue
        operator = getproperty(model.operators, operator_name)
        number_complex_values += prod(size(operator.coefficients)[1:(end - 1)])
    end
    bytes_per_kpoint = sizeof(T) * number_complex_values
    memory_limited = max(1, _DEFAULT_INTERPOLATION_MEMORY ÷ max(1, bytes_per_kpoint))
    return min(number_kpoints, memory_limited, _MAXIMUM_INTERPOLATION_BATCH)
end

function _fourier_phase_block!(
        phase::AbstractMatrix{T},
        real_space::RealSpaceDomain,
        kpoints::AbstractVector,
    ) where {T <: Complex}
    size(phase) == (length(real_space), length(kpoints)) ||
        throw(DimensionMismatch("phase block has the wrong shape"))
    for (batch_index, kpoint) in enumerate(kpoints)
        length(kpoint) == 3 || throw(ArgumentError("each k point must have length three"))
        for (vector_index, vector) in enumerate(real_space.vectors)
            phase[vector_index, batch_index] = cis(2π * dot(kpoint, vector))
        end
    end
    return phase
end

function _fourier_phase_block(
        real_space::RealSpaceDomain,
        kpoints::AbstractVector,
        kpoint_range::AbstractUnitRange,
        ::Type{T},
    ) where {T <: Complex}
    phase = Matrix{T}(undef, length(real_space), length(kpoint_range))
    return _fourier_phase_block!(phase, real_space, view(kpoints, kpoint_range))
end

function _evaluate_real_space_operator!(
        destination::AbstractArray,
        operator::RealSpaceOperator,
        phase::AbstractMatrix,
    )
    coefficients = operator.coefficients
    number_vectors = size(coefficients, ndims(coefficients))
    size(phase, 1) == number_vectors ||
        throw(DimensionMismatch("operator and phase-block domains differ"))
    size(destination)[1:(end - 1)] == size(coefficients)[1:(end - 1)] ||
        throw(DimensionMismatch("operator and destination value shapes differ"))
    size(destination, ndims(destination)) == size(phase, 2) ||
        throw(DimensionMismatch("destination and phase-block batch sizes differ"))

    packed_coefficients = reshape(coefficients, :, number_vectors)
    packed_destination = reshape(destination, :, size(phase, 2))
    mul!(packed_destination, packed_coefficients, phase)
    return destination
end

function _evaluate_real_space_operator(
        operator::RealSpaceOperator,
        real_space::RealSpaceDomain,
        kpoints::AbstractVector,
        kpoint_range::AbstractUnitRange,
    )
    coefficients = operator.coefficients
    T = eltype(coefficients)
    phase = _fourier_phase_block(real_space, kpoints, kpoint_range, T)
    values = Array{T}(
        undef, size(coefficients)[1:(end - 1)]..., length(kpoint_range)
    )
    return _evaluate_real_space_operator!(values, operator, phase)
end

function _evaluate_real_space_operator_gradient!(
        destination::AbstractArray,
        operator::RealSpaceOperator,
        real_space::RealSpaceDomain,
        phase::AbstractMatrix,
        derivative_phase::AbstractMatrix,
        derivative_values::AbstractMatrix,
    )
    number_wannier = n_wannier(operator)
    number_vectors = n_Rvectors(operator)
    number_kpoints = size(phase, 2)
    value_shape = size(operator.coefficients)[1:(end - 1)]
    size(destination) == (value_shape..., 3, number_kpoints) || throw(
        DimensionMismatch("operator-gradient destination has the wrong shape"),
    )
    size(derivative_phase) == size(phase) ||
        throw(DimensionMismatch("derivative and ordinary phase blocks differ"))
    number_values = prod(value_shape)
    size(derivative_values, 1) >= number_values &&
        size(derivative_values, 2) == number_kpoints ||
        throw(DimensionMismatch("operator-gradient buffer has the wrong shape"))

    packed_coefficients = reshape(operator.coefficients, number_values, number_vectors)
    active_derivative_values = view(derivative_values, 1:number_values, :)
    for cartesian_component in 1:3
        for kpoint_index in 1:number_kpoints, vector_index in 1:number_vectors
            derivative_phase[vector_index, kpoint_index] =
                im * real_space.cartesian_vectors[vector_index][cartesian_component] *
                phase[vector_index, kpoint_index]
        end
        mul!(
            active_derivative_values,
            packed_coefficients,
            derivative_phase,
        )
        component_destination = selectdim(
            destination, ndims(destination) - 1, cartesian_component
        )
        copyto!(
            component_destination, reshape(
                active_derivative_values, value_shape..., number_kpoints
            )
        )
    end
    return destination
end

function _evaluate_hamiltonian_gradient!(
        destination::AbstractArray{<:Complex, 4},
        operator::RealSpaceOperator,
        real_space::RealSpaceDomain,
        phase::AbstractMatrix,
        derivative_phase::AbstractMatrix,
        derivative_values::AbstractMatrix,
    )
    component_shape(operator) == () ||
        throw(ArgumentError("Hamiltonian derivatives require a scalar operator"))
    return _evaluate_real_space_operator_gradient!(
        destination,
        operator,
        real_space,
        phase,
        derivative_phase,
        derivative_values,
    )
end
