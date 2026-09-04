function _validate_bloch_operator(
        name::Symbol,
        operator::BlochOperator,
        number_bands::Integer,
        number_kpoints::Integer,
    )
    values = operator.values
    if ndims(values) == 2
        size(values) == (number_bands, number_kpoints) || throw(
            ArgumentError(
                "diagonal Bloch operator :$name must have size " *
                    "($number_bands, $number_kpoints), got $(size(values))",
            ),
        )
        isempty(component_shape(operator.law)) || throw(
            ArgumentError("a diagonal Bloch operator can only have a scalar transformation law"),
        )
        return nothing
    end

    size(values, 1) == number_bands && size(values, 2) == number_bands || throw(
        ArgumentError(
            "Bloch operator :$name must have two leading band axes of length $number_bands",
        ),
    )
    size(values, ndims(values)) == number_kpoints ||
        throw(ArgumentError("the final axis of Bloch operator :$name must index k points"))
    actual_shape = Tuple(size(values)[3:(end - 1)])
    expected_shape = component_shape(operator.law)
    actual_shape == expected_shape || throw(
        ArgumentError(
            "Bloch operator :$name has component shape $actual_shape; " *
                "its transformation law requires $expected_shape",
        ),
    )
    return nothing
end

function _transform_bloch_operator(operator::BlochOperator, gauges::AbstractArray{<:Complex, 3})
    values = operator.values
    number_bands, number_wannier, number_kpoints = size(gauges)
    T = promote_type(eltype(values), eltype(gauges))

    if ndims(values) == 2
        transformed = Array{T}(undef, number_wannier, number_wannier, number_kpoints)
        for kpoint_index in 1:number_kpoints
            gauge = view(gauges, :, :, kpoint_index)
            diagonal = Diagonal(view(values, :, kpoint_index))
            rotated = gauge' * diagonal * gauge
            if operator.hermitian
                transformed[:, :, kpoint_index] .= Hermitian(rotated)
            else
                transformed[:, :, kpoint_index] .= rotated
            end
        end
        return transformed
    end

    components = component_shape(operator.law)
    number_components = prod(components; init = 1)
    source = reshape(values, number_bands, number_bands, number_components, number_kpoints)
    transformed = Array{T}(
        undef, number_wannier, number_wannier, number_components, number_kpoints
    )
    for kpoint_index in 1:number_kpoints
        gauge = view(gauges, :, :, kpoint_index)
        for component_index in 1:number_components
            matrix = view(source, :, :, component_index, kpoint_index)
            rotated = gauge' * matrix * gauge
            if operator.hermitian
                transformed[:, :, component_index, kpoint_index] .= Hermitian(rotated)
            else
                transformed[:, :, component_index, kpoint_index] .= rotated
            end
        end
    end
    return reshape(
        transformed, number_wannier, number_wannier, components..., number_kpoints
    )
end

function _quotient_fourier_phase(
        kpoints::AbstractVector,
        vectors::AbstractVector,
        ::Type{T},
    ) where {T <: Complex}
    number_kpoints = length(kpoints)
    phase = Matrix{T}(undef, number_kpoints, length(vectors))
    normalization = inv(number_kpoints)
    for (vector_index, vector) in enumerate(vectors)
        for (kpoint_index, kpoint) in enumerate(kpoints)
            phase[kpoint_index, vector_index] =
                cis(-2π * dot(kpoint, vector)) * normalization
        end
    end
    return phase
end

function _quotient_fourier_coefficients(values::AbstractArray, phase::AbstractMatrix)
    number_kpoints = size(values, ndims(values))
    size(phase, 1) == number_kpoints ||
        throw(DimensionMismatch("the Fourier phase and operator k-point axes differ"))
    packed_values = reshape(values, :, number_kpoints)
    T = promote_type(eltype(values), eltype(phase))
    packed_coefficients = Matrix{T}(undef, size(packed_values, 1), size(phase, 2))
    mul!(packed_coefficients, packed_values, phase)
    return reshape(packed_coefficients, size(values)[1:(end - 1)]..., size(phase, 2))
end

function _canonical_vector_order(real_space::RealSpaceDomain, vectors::AbstractVector)
    input_index = Dict(Vec3{Int}(vector) => index for (index, vector) in enumerate(vectors))
    return map(real_space.vectors) do vector
        get(input_index, vector) do
            error("canonical real-space vector is absent from its source domain")
        end
    end
end

function _reorder_vector_axis(coefficients::AbstractArray, order::AbstractVector{<:Integer})
    indices = (ntuple(_ -> Colon(), ndims(coefficients) - 1)..., order)
    return coefficients[indices...]
end

function _pack_real_space_operators(
        lattice,
        representative_vectors,
        operator_descriptions::NamedTuple,
        selected_coefficients,
    )
    domain = RealSpaceDomain(lattice, representative_vectors)
    order = _canonical_vector_order(domain, representative_vectors)
    real_space_operators = map(
        selected_coefficients, values(operator_descriptions)
    ) do coefficients, description
        sorted_coefficients = _reorder_vector_axis(coefficients, order)
        RealSpaceOperator(
            sorted_coefficients,
            description.law,
            domain;
            hermitian = description.hermitian,
        )
    end
    return domain, NamedTuple{keys(operator_descriptions)}(real_space_operators)
end

function _interpolation_crystal(lattice, atom_positions, atom_labels)
    length(atom_positions) == length(atom_labels) ||
        throw(DimensionMismatch("atom-position and atom-label counts differ"))
    T = float(eltype(lattice))
    return InterpolationCrystal{T}(
        Mat3{T}(lattice), Vector{Vec3{T}}(atom_positions), String.(atom_labels)
    )
end

function _wannier_basis(fractional_centers, number_wannier::Integer, ::Type{T}) where {T}
    length(fractional_centers) == number_wannier || throw(
        DimensionMismatch(
            "expected $number_wannier fractional Wannier centers, " *
                "got $(length(fractional_centers))",
        ),
    )
    return WannierBasis{T}(Vector{Vec3{T}}(fractional_centers))
end

function InterpolationModel(
        model::Model;
        operators::NamedTuple = (;),
        real_space::RealSpaceScheme = MinimumDistance(),
        symmetry = nothing,
    )
    isnothing(symmetry) || symmetry isa WannierSymmetry ||
        throw(ArgumentError("symmetry must be a WannierSymmetry or nothing"))
    haskey(operators, :hamiltonian) &&
        throw(ArgumentError(":hamiltonian is constructed from model.eigenvalues"))

    hamiltonian = BlochOperator(
        model.eigenvalues; law = Scalar(time_reversal = Even()), hermitian = true
    )
    bloch_operators = merge((; hamiltonian), operators)
    number_bands = n_bands(model)
    number_kpoints = n_kpoints(model)
    for (name, operator) in pairs(bloch_operators)
        operator isa BlochOperator ||
            throw(ArgumentError("operator :$name must be a BlochOperator"))
        _validate_bloch_operator(name, operator, number_bands, number_kpoints)
    end

    fractional_centers = if isnothing(symmetry)
        inverse_lattice = inv(model.lattice)
        map(center(model)) do cartesian_center
            Vec3(inverse_lattice * cartesian_center)
        end
    else
        n_wannier(symmetry) == n_wannier(model) || throw(
            DimensionMismatch("model and Wannier symmetry have different basis sizes"),
        )
        symmetry.centers
    end
    selection = _real_space_selection(
        real_space, model.lattice, kgrid_size(model), fractional_centers
    )

    transformed_operators = map(values(bloch_operators)) do operator
        _transform_bloch_operator(operator, model.gauges)
    end
    coefficient_type = promote_type(map(eltype, transformed_operators)...)
    phase = _quotient_fourier_phase(
        kpoints(model), _quotient_vectors(selection), coefficient_type
    )
    selected_coefficients = map(transformed_operators) do values
        quotient_coefficients = _quotient_fourier_coefficients(values, phase)
        _apply_real_space_selection(selection, quotient_coefficients)
    end

    representative_vectors = _representative_vectors(selection)
    if !isnothing(symmetry)
        representative_vectors, selected_coefficients = _close_real_space_operators(
            model.lattice,
            representative_vectors,
            bloch_operators,
            selected_coefficients,
            symmetry,
        )
    end

    domain, operator_tuple = _pack_real_space_operators(
        model.lattice, representative_vectors, bloch_operators, selected_coefficients
    )
    crystal = _interpolation_crystal(
        model.lattice, model.atom_positions, model.atom_labels
    )
    basis = _wannier_basis(
        fractional_centers, n_wannier(model), eltype(crystal.lattice)
    )
    return InterpolationModel(crystal, basis, domain, operator_tuple, symmetry)
end
