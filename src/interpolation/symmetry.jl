function _real_space_dictionary(
        coefficients::AbstractArray,
        vectors::AbstractVector,
    )
    ndims(coefficients) >= 3 ||
        throw(ArgumentError("real-space operator coefficients need two matrix axes"))
    size(coefficients, ndims(coefficients)) == length(vectors) ||
        throw(DimensionMismatch("coefficient and real-space-vector counts differ"))
    number_wannier = size(coefficients, 1)
    number_components = prod(size(coefficients)[3:(end - 1)]; init = 1)
    packed = reshape(
        coefficients, number_wannier, number_wannier, number_components, length(vectors)
    )
    return Dict(
        Vec3{Int}(vector) => Array(view(packed, :, :, :, index))
            for (index, vector) in enumerate(vectors)
    )
end

function _zero_real_space_matrix(realization::AbstractDict)
    first_matrix = first(values(realization))
    return zeros(eltype(first_matrix), size(first_matrix))
end

function _hermitian_real_space_projection(realization::AbstractDict)
    vectors = union(collect(keys(realization)), [-vector for vector in keys(realization)])
    zero_matrix = _zero_real_space_matrix(realization)
    projected = empty(realization)
    for vector in vectors
        forward = get(realization, vector, zero_matrix)
        backward = get(realization, -vector, zero_matrix)
        value = similar(forward)
        for component_index in axes(forward, 3)
            view(value, :, :, component_index) .=
                (
                view(forward, :, :, component_index) +
                    view(backward, :, :, component_index)'
            ) / 2
        end
        projected[vector] = value
    end
    return projected
end

_time_reversal_sign(::Even) = 1
_time_reversal_sign(::Odd) = -1

_time_reversal_parity(law::OperatorLaw) = law.time_reversal

_supports_homogeneous_symmetry_closure(::OperatorLaw) = true
_supports_homogeneous_symmetry_closure(::_HamiltonianPositionLaw) = false
_supports_homogeneous_symmetry_closure(::_PositionHamiltonianPositionLaw) = false

function _shift_position_center!(
        realization::AbstractDict,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
        sign::Integer,
    )
    origin = Vec3(0, 0, 0)
    origin_matrix = get!(realization, origin) do
        _zero_real_space_matrix(realization)
    end
    cartesian_centers = map(center -> lattice * center, symmetry.centers)
    for band_index in eachindex(cartesian_centers), component_index in 1:3
        origin_matrix[band_index, band_index, component_index] +=
            sign * cartesian_centers[band_index][component_index]
    end
    return realization
end

function _close_centered_position(
        vectors::AbstractVector,
        coefficients::AbstractArray,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    realization = _real_space_dictionary(coefficients, vectors)
    _shift_position_center!(realization, symmetry, lattice, -1)
    projected = _close_homogeneous_realization(
        realization,
        PolarVector(time_reversal = Even()),
        true,
        symmetry,
        lattice,
    )
    _shift_position_center!(projected, symmetry, lattice, 1)
    return projected
end

function _close_homogeneous_realization(
        realization::AbstractDict,
        law::OperatorLaw,
        hermitian::Bool,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    hermitian && (realization = _hermitian_real_space_projection(realization))
    projected = _symmetry_project_real_space(realization, law, symmetry, lattice)
    hermitian && (projected = _hermitian_real_space_projection(projected))
    return projected
end

function _cartesian_rotation(operation::SymOp, lattice::AbstractMatrix)
    rotation = lattice * operation.W * inv(lattice)
    return real.(rotation)
end

function _component_representation(
        ::Scalar,
        ::SymOp,
        ::AbstractMatrix,
    )
    return ones(1, 1)
end

function _component_representation(
        ::PolarVector,
        operation::SymOp,
        lattice::AbstractMatrix,
    )
    return _cartesian_rotation(operation, lattice)
end

function _component_representation(
        ::AxialVector,
        operation::SymOp,
        lattice::AbstractMatrix,
    )
    rotation = _cartesian_rotation(operation, lattice)
    return sign(det(rotation)) .* rotation
end

function _component_representation(
        law::CartesianTensor,
        operation::SymOp,
        lattice::AbstractMatrix,
    )
    law.rank == 0 && return ones(1, 1)
    rotation = _cartesian_rotation(operation, lattice)
    component_indices = CartesianIndices(component_shape(law))
    representation = Matrix{eltype(rotation)}(
        undef, length(component_indices), length(component_indices)
    )
    axial_factor = law.axial ? sign(det(rotation)) : one(eltype(rotation))
    for (output_linear, output_index) in enumerate(component_indices)
        for (input_linear, input_index) in enumerate(component_indices)
            value = axial_factor
            for tensor_axis in 1:law.rank
                value *= rotation[output_index[tensor_axis], input_index[tensor_axis]]
            end
            representation[output_linear, input_linear] = value
        end
    end
    return representation
end

function _symmetry_project_real_space(
        realization::AbstractDict,
        law::OperatorLaw,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    number_wannier = n_wannier(symmetry)
    number_components = prod(component_shape(law); init = 1)
    number_symmetries = length(symmetry.symops)
    projected = empty(realization)

    for (symmetry_index, operation) in enumerate(symmetry.symops)
        representation = symmetry.orbital_reps[symmetry_index].D
        translations = symmetry.translations[symmetry_index]
        targets = [findall(!iszero, view(representation, :, source)) for source in 1:number_wannier]
        component_representation = _component_representation(law, operation, lattice)
        size(component_representation) == (number_components, number_components) ||
            error("operator law returned a component representation of the wrong size")
        parity = operation.time_reversal ?
            _time_reversal_sign(_time_reversal_parity(law)) : 1

        for (vector, matrix) in realization
            transformed_vector = operation.W * vector
            for column in 1:number_wannier, row in 1:number_wannier
                output_vector = Vec3{Int}(
                    transformed_vector + translations[column] - translations[row]
                )
                output_matrix = get!(projected, output_vector) do
                    zeros(
                        eltype(matrix), number_wannier, number_wannier, number_components
                    )
                end
                for input_component in 1:number_components
                    value = matrix[row, column, input_component]
                    operation.time_reversal && (value = conj(value))
                    iszero(value) && continue
                    for output_component in 1:number_components
                        component_factor = component_representation[
                            output_component, input_component,
                        ]
                        iszero(component_factor) && continue
                        for output_column in targets[column], output_row in targets[row]
                            output_matrix[output_row, output_column, output_component] +=
                                parity * component_factor *
                                representation[output_row, row] * value *
                                conj(representation[output_column, column]) /
                                number_symmetries
                        end
                    end
                end
            end
        end
    end
    return projected
end

function _close_real_space_operator(
        vectors::AbstractVector,
        description,
        coefficients::AbstractArray,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    description.law isa _BerryConnectionLaw && return _close_centered_position(
        vectors, coefficients, symmetry, lattice
    )
    _supports_homogeneous_symmetry_closure(description.law) || throw(
        ArgumentError(
            "symmetry closure of $(nameof(typeof(description.law))) requires " *
                "its dedicated inhomogeneous transformation law and is not " *
                "implemented yet",
        ),
    )
    size(coefficients, 1) == n_wannier(symmetry) || throw(
        DimensionMismatch("operator and Wannier symmetry have different basis sizes"),
    )

    realization = _real_space_dictionary(coefficients, vectors)
    return _close_homogeneous_realization(
        realization,
        description.law,
        description.hermitian,
        symmetry,
        lattice,
    )
end

function _realization_vector_union(realizations::AbstractDict...; hermitian = false)
    vectors = Vec3{Int}[]
    for realization in realizations
        append!(vectors, keys(realization))
        hermitian && append!(vectors, (-vector for vector in keys(realization)))
    end
    return unique!(vectors)
end

function _center_hamiltonian_position(
        hamiltonian::AbstractDict,
        hamiltonian_position::AbstractDict,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    cartesian_centers = map(center -> lattice * center, symmetry.centers)
    zero_hamiltonian = _zero_real_space_matrix(hamiltonian)
    zero_position = _zero_real_space_matrix(hamiltonian_position)
    centered = empty(hamiltonian_position)
    for vector in _realization_vector_union(hamiltonian, hamiltonian_position)
        hamiltonian_matrix = get(hamiltonian, vector, zero_hamiltonian)
        position_matrix = get(hamiltonian_position, vector, zero_position)
        value = copy(position_matrix)
        for column in axes(value, 2), row in axes(value, 1), component_index in 1:3
            value[row, column, component_index] -=
                hamiltonian_matrix[row, column, 1] *
                cartesian_centers[column][component_index]
        end
        centered[vector] = value
    end
    return centered
end

function _restore_hamiltonian_position(
        hamiltonian::AbstractDict,
        centered_position::AbstractDict,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    cartesian_centers = map(center -> lattice * center, symmetry.centers)
    zero_hamiltonian = _zero_real_space_matrix(hamiltonian)
    zero_position = _zero_real_space_matrix(centered_position)
    restored = empty(centered_position)
    for vector in _realization_vector_union(hamiltonian, centered_position)
        hamiltonian_matrix = get(hamiltonian, vector, zero_hamiltonian)
        centered_matrix = get(centered_position, vector, zero_position)
        value = copy(centered_matrix)
        for column in axes(value, 2), row in axes(value, 1), component_index in 1:3
            value[row, column, component_index] +=
                hamiltonian_matrix[row, column, 1] *
                cartesian_centers[column][component_index]
        end
        restored[vector] = value
    end
    return restored
end

_tensor_component(first_component::Integer, second_component::Integer) =
    first_component + 3(second_component - 1)

function _center_position_hamiltonian_position(
        hamiltonian::AbstractDict,
        hamiltonian_position::AbstractDict,
        position_hamiltonian_position::AbstractDict,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    cartesian_centers = map(center -> lattice * center, symmetry.centers)
    zero_hamiltonian = _zero_real_space_matrix(hamiltonian)
    zero_position = _zero_real_space_matrix(hamiltonian_position)
    zero_second_moment = _zero_real_space_matrix(position_hamiltonian_position)
    centered = empty(position_hamiltonian_position)
    vectors = _realization_vector_union(
        hamiltonian,
        hamiltonian_position,
        position_hamiltonian_position;
        hermitian = true,
    )
    for vector in vectors
        hamiltonian_matrix = get(hamiltonian, vector, zero_hamiltonian)
        position_matrix = get(hamiltonian_position, vector, zero_position)
        adjoint_position_matrix = get(
            hamiltonian_position, -vector, zero_position
        )
        second_moment_matrix = get(
            position_hamiltonian_position, vector, zero_second_moment
        )
        value = copy(second_moment_matrix)
        for column in axes(value, 2), row in axes(value, 1)
            row_center = cartesian_centers[row]
            column_center = cartesian_centers[column]
            for second_component in 1:3, first_component in 1:3
                component_index = _tensor_component(
                    first_component, second_component
                )
                left_position = conj(
                    adjoint_position_matrix[column, row, first_component]
                )
                value[row, column, component_index] -=
                    left_position * column_center[second_component] +
                    row_center[first_component] *
                    position_matrix[row, column, second_component] -
                    row_center[first_component] *
                    hamiltonian_matrix[row, column, 1] *
                    column_center[second_component]
            end
        end
        centered[vector] = value
    end
    return centered
end

function _restore_position_hamiltonian_position(
        hamiltonian::AbstractDict,
        centered_hamiltonian_position::AbstractDict,
        centered_second_moment::AbstractDict,
        symmetry::WannierSymmetry,
        lattice::AbstractMatrix,
    )
    cartesian_centers = map(center -> lattice * center, symmetry.centers)
    zero_hamiltonian = _zero_real_space_matrix(hamiltonian)
    zero_position = _zero_real_space_matrix(centered_hamiltonian_position)
    zero_second_moment = _zero_real_space_matrix(centered_second_moment)
    restored = empty(centered_second_moment)
    vectors = _realization_vector_union(
        hamiltonian,
        centered_hamiltonian_position,
        centered_second_moment;
        hermitian = true,
    )
    for vector in vectors
        hamiltonian_matrix = get(hamiltonian, vector, zero_hamiltonian)
        position_matrix = get(
            centered_hamiltonian_position, vector, zero_position
        )
        adjoint_position_matrix = get(
            centered_hamiltonian_position, -vector, zero_position
        )
        second_moment_matrix = get(
            centered_second_moment, vector, zero_second_moment
        )
        value = copy(second_moment_matrix)
        for column in axes(value, 2), row in axes(value, 1)
            row_center = cartesian_centers[row]
            column_center = cartesian_centers[column]
            for second_component in 1:3, first_component in 1:3
                component_index = _tensor_component(
                    first_component, second_component
                )
                left_position = conj(
                    adjoint_position_matrix[column, row, first_component]
                )
                value[row, column, component_index] +=
                    left_position * column_center[second_component] +
                    row_center[first_component] *
                    position_matrix[row, column, second_component] +
                    row_center[first_component] *
                    hamiltonian_matrix[row, column, 1] *
                    column_center[second_component]
            end
        end
        restored[vector] = value
    end
    return restored
end

function _second_moment_hermitian_projection(realization::AbstractDict)
    vectors = _realization_vector_union(realization; hermitian = true)
    zero_matrix = _zero_real_space_matrix(realization)
    projected = empty(realization)
    for vector in vectors
        forward = get(realization, vector, zero_matrix)
        backward = get(realization, -vector, zero_matrix)
        value = similar(forward)
        for column in axes(value, 2), row in axes(value, 1)
            for second_component in 1:3, first_component in 1:3
                component_index = _tensor_component(
                    first_component, second_component
                )
                transposed_component = _tensor_component(
                    second_component, first_component
                )
                value[row, column, component_index] = (
                    forward[row, column, component_index] +
                        conj(backward[column, row, transposed_component])
                ) / 2
            end
        end
        projected[vector] = value
    end
    return projected
end

function _pack_real_space_dictionary(
        realization::AbstractDict,
        vectors::AbstractVector,
        law::OperatorLaw,
    )
    first_matrix = first(values(realization))
    number_wannier = size(first_matrix, 1)
    packed = zeros(
        eltype(first_matrix),
        number_wannier,
        number_wannier,
        size(first_matrix, 3),
        length(vectors),
    )
    for (vector_index, vector) in enumerate(vectors)
        haskey(realization, vector) || continue
        view(packed, :, :, :, vector_index) .= realization[vector]
    end
    return reshape(
        packed, number_wannier, number_wannier, component_shape(law)..., length(vectors)
    )
end

function _close_real_space_operators(
        lattice::AbstractMatrix,
        vectors::AbstractVector,
        descriptions::NamedTuple,
        coefficients,
        symmetry::WannierSymmetry,
    )
    names = keys(descriptions)
    coefficient_tuple = NamedTuple{names}(coefficients)
    raw_realizations = NamedTuple{names}(
        map(values -> _real_space_dictionary(values, vectors), coefficient_tuple)
    )
    closed_by_name = Dict{Symbol, Any}()
    for name in names
        description = getproperty(descriptions, name)
        description.law isa _HamiltonianPositionLaw && continue
        description.law isa _PositionHamiltonianPositionLaw && continue
        closed_by_name[name] = _close_real_space_operator(
            vectors,
            description,
            getproperty(coefficient_tuple, name),
            symmetry,
            lattice,
        )
    end

    if haskey(descriptions, :hamiltonian_position)
        haskey(descriptions, :hamiltonian) || throw(
            ArgumentError("HamiltonianPosition symmetry closure requires :hamiltonian"),
        )
        centered_position = _center_hamiltonian_position(
            raw_realizations.hamiltonian,
            raw_realizations.hamiltonian_position,
            symmetry,
            lattice,
        )
        centered_position = _close_homogeneous_realization(
            centered_position,
            PolarVector(time_reversal = Even()),
            false,
            symmetry,
            lattice,
        )
        closed_by_name[:hamiltonian_position] = _restore_hamiltonian_position(
            closed_by_name[:hamiltonian],
            centered_position,
            symmetry,
            lattice,
        )
    end

    if haskey(descriptions, :position_hamiltonian_position)
        haskey(descriptions, :hamiltonian_position) || throw(
            ArgumentError(
                "PositionHamiltonianPosition symmetry closure requires " *
                    ":hamiltonian_position",
            ),
        )
        centered_second_moment = _center_position_hamiltonian_position(
            raw_realizations.hamiltonian,
            raw_realizations.hamiltonian_position,
            raw_realizations.position_hamiltonian_position,
            symmetry,
            lattice,
        )
        centered_second_moment = _second_moment_hermitian_projection(
            centered_second_moment
        )
        centered_second_moment = _symmetry_project_real_space(
            centered_second_moment,
            CartesianTensor(2; time_reversal = Even()),
            symmetry,
            lattice,
        )
        centered_second_moment = _second_moment_hermitian_projection(
            centered_second_moment
        )
        centered_position = _center_hamiltonian_position(
            closed_by_name[:hamiltonian],
            closed_by_name[:hamiltonian_position],
            symmetry,
            lattice,
        )
        closed_by_name[:position_hamiltonian_position] =
            _restore_position_hamiltonian_position(
            closed_by_name[:hamiltonian],
            centered_position,
            centered_second_moment,
            symmetry,
            lattice,
        )
    end
    realizations = map(name -> closed_by_name[name], names)
    closed_vectors = sort!(
        unique!(vcat((collect(keys(realization)) for realization in realizations)...));
        by = Tuple,
    )
    packed = map(realizations, descriptions) do realization, description
        _pack_real_space_dictionary(realization, closed_vectors, description.law)
    end
    return closed_vectors, packed
end

"""
Independent query-point Reynolds projection for lattice-periodic operators.
This is a reference implementation for validating one-shot real-space closure;
production interpolation evaluates the already-closed coefficients directly.
"""
function _project_operator_at_kpoint(
        operator::RealSpaceOperator,
        real_space::RealSpaceDomain,
        lattice::AbstractMatrix,
        symmetry::WannierSymmetry,
        kpoint::AbstractVector,
    )
    number_wannier = n_wannier(operator)
    number_wannier == n_wannier(symmetry) || throw(
        DimensionMismatch("operator and Wannier symmetry have different basis sizes"),
    )
    number_components = prod(component_shape(operator.law); init = 1)
    result = zeros(
        eltype(operator.coefficients), number_wannier, number_wannier, number_components
    )
    for (symmetry_index, operation) in enumerate(symmetry.symops)
        inverse_operation = symmetry.symops[operation.isym_inv]
        source_kpoint = rotate_kpoint(kpoint, inverse_operation)
        source_logical = _evaluate_real_space_operator(
            operator, real_space, [source_kpoint], Base.OneTo(1)
        )
        source = reshape(source_logical, number_wannier, number_wannier, number_components)
        transport = _expansion_matrix(symmetry_index, source_kpoint, symmetry)
        component_representation = _component_representation(
            operator.law, operation, lattice
        )
        parity = operation.time_reversal ?
            _time_reversal_sign(_time_reversal_parity(operator.law)) : 1
        for input_component in 1:number_components
            transformed = transport' * view(source, :, :, input_component) * transport
            operation.time_reversal && (transformed = conj.(transformed))
            for output_component in 1:number_components
                view(result, :, :, output_component) .+=
                    parity * component_representation[output_component, input_component] .*
                    transformed ./ length(symmetry.symops)
            end
        end
    end
    return reshape(result, number_wannier, number_wannier, component_shape(operator.law)...)
end

function _symmetry_covariance_residual(
        model::InterpolationModel,
        operator_name::Symbol,
        kpoints::AbstractVector,
    )
    isnothing(model.symmetry) &&
        throw(ArgumentError("symmetry diagnostics require a symmetry-closed model"))
    haskey(model.operators, operator_name) ||
        throw(ArgumentError("interpolation model has no operator :$operator_name"))
    operator = getproperty(model.operators, operator_name)
    number_wannier = n_wannier(operator)
    number_components = prod(component_shape(operator.law); init = 1)
    maximum_residual = 0.0
    worst_kpoint = 0
    worst_symmetry = 0

    for (kpoint_index, kpoint) in enumerate(kpoints)
        source_logical = _evaluate_real_space_operator(
            operator, model.real_space, [kpoint], Base.OneTo(1)
        )
        source = reshape(source_logical, number_wannier, number_wannier, number_components)
        for (symmetry_index, operation) in enumerate(model.symmetry.symops)
            target_kpoint = rotate_kpoint(kpoint, operation)
            target_logical = _evaluate_real_space_operator(
                operator, model.real_space, [target_kpoint], Base.OneTo(1)
            )
            target = reshape(
                target_logical, number_wannier, number_wannier, number_components
            )
            transport = _expansion_matrix(symmetry_index, kpoint, model.symmetry)
            component_representation = _component_representation(
                operator.law, operation, model.crystal.lattice
            )
            parity = operation.time_reversal ?
                _time_reversal_sign(_time_reversal_parity(operator.law)) : 1
            expected = zeros(eltype(source), size(source))
            for input_component in 1:number_components
                transformed = transport' * view(source, :, :, input_component) * transport
                operation.time_reversal && (transformed = conj.(transformed))
                for output_component in 1:number_components
                    view(expected, :, :, output_component) .+=
                        parity *
                        component_representation[output_component, input_component] .*
                        transformed
                end
            end
            residual = norm(target - expected)
            if residual > maximum_residual
                maximum_residual = residual
                worst_kpoint = kpoint_index
                worst_symmetry = symmetry_index
            end
        end
    end
    return (;
        maximum = maximum_residual,
        kpoint_index = worst_kpoint,
        symmetry_index = worst_symmetry,
    )
end
