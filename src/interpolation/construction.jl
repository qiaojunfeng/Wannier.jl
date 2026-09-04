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

function _reduce_real_space_coefficients(reducer, coefficients::AbstractArray)
    number_dimensions = ndims(coefficients)
    number_wannier = size(coefficients, 1)
    number_vectors = size(coefficients, number_dimensions)
    components = Tuple(size(coefficients)[3:(end - 1)])
    number_components = prod(components; init = 1)
    packed = reshape(
        coefficients, number_wannier, number_wannier, number_components, number_vectors
    )

    number_reduced_vectors = length(reducer.Rvectors)
    reduced = Array{eltype(coefficients)}(
        undef, number_wannier, number_wannier, number_components, number_reduced_vectors
    )
    for component_index in 1:number_components
        reduced[:, :, component_index, :] .=
            reducer(view(packed, :, :, component_index, :))
    end
    return reshape(
        reduced,
        number_wannier,
        number_wannier,
        components...,
        number_reduced_vectors,
    )
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

function _legacy_real_space(
        model::Model,
        scheme::WignerSeitz,
        fractional_centers,
    )
    return WignerSeitzRspace(
        model.lattice,
        kgrid_size(model);
        atol = scheme.atol,
        max_cell = scheme.max_cell,
    )
end

function _legacy_real_space(
        model::Model,
        scheme::MinimumDistance,
        fractional_centers,
    )
    wigner_seitz = WignerSeitzRspace(
        model.lattice,
        kgrid_size(model);
        atol = scheme.atol,
        max_cell = scheme.max_cell,
    )
    return MDRSRspace(
        wigner_seitz,
        kgrid_size(model),
        fractional_centers;
        atol = scheme.atol,
        max_cell = scheme.max_cell,
    )
end

function InterpolationModel(
        model::Model;
        operators::NamedTuple = (;),
        real_space::RealSpaceScheme = MinimumDistance(),
        symmetry = nothing,
    )
    isnothing(symmetry) || throw(
        ArgumentError(
            "symmetry-closed real-space construction is not implemented in this first redesign slice",
        ),
    )
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

    inverse_lattice = inv(model.lattice)
    fractional_centers = map(center(model)) do cartesian_center
        Vec3(inverse_lattice * cartesian_center)
    end
    selected_real_space = _legacy_real_space(model, real_space, fractional_centers)
    reducer = RvectorReducer(selected_real_space)

    transformed_operators = map(values(bloch_operators)) do operator
        _transform_bloch_operator(operator, model.gauges)
    end
    coefficient_type = promote_type(map(eltype, transformed_operators)...)
    phase = _quotient_fourier_phase(
        kpoints(model), selected_real_space.Rvectors, coefficient_type
    )
    reduced_coefficients = map(transformed_operators) do values
        quotient_coefficients = _quotient_fourier_coefficients(values, phase)
        _reduce_real_space_coefficients(reducer, quotient_coefficients)
    end

    domain = RealSpaceDomain(model.lattice, reducer.Rvectors)
    order = _canonical_vector_order(domain, reducer.Rvectors)
    real_space_operators = map(reduced_coefficients, values(bloch_operators)) do coefficients, input
        sorted_coefficients = _reorder_vector_axis(coefficients, order)
        RealSpaceOperator(sorted_coefficients, input.law, domain)
    end
    operator_names = keys(bloch_operators)
    operator_tuple = NamedTuple{operator_names}(real_space_operators)

    T = eltype(model.lattice)
    crystal = InterpolationCrystal{T}(
        model.lattice, copy(model.atom_positions), copy(model.atom_labels)
    )
    basis = WannierBasis{T}(Vector{Vec3{T}}(fractional_centers))
    return InterpolationModel(crystal, basis, domain, operator_tuple, symmetry)
end
