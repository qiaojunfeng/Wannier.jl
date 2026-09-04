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
    size(coefficients, 1) == n_wannier(symmetry) || throw(
        DimensionMismatch("operator and Wannier symmetry have different basis sizes"),
    )

    realization = _real_space_dictionary(coefficients, vectors)
    description.hermitian && (realization = _hermitian_real_space_projection(realization))
    projected = _symmetry_project_real_space(
        realization, description.law, symmetry, lattice
    )
    description.hermitian && (projected = _hermitian_real_space_projection(projected))
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
    realizations = map(descriptions, coefficients) do description, values
        _close_real_space_operator(vectors, description, values, symmetry, lattice)
    end
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
