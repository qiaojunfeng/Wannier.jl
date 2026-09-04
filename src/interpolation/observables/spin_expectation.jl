export SpinExpectation

"""
    SpinExpectation()
    SpinExpectation(direction; truncate = false)

Request diagonal band-basis expectation values of the primitive `:spin`
operator. Without a direction, `result.spin_expectation` has layout
`3 × n_wannier × n_kpoints`. With a Cartesian direction, it is projected onto
that normalized axis and has layout `n_wannier × n_kpoints`.

Values have the units of the primitive operator and are dimensionless for a
Wannier90 `spn` file. At exact degeneracies, individual diagonal expectations
depend on the eigenvector basis selected within the degenerate subspace. Set
`truncate = true` to clamp numerical values to `[-1, 1]`.
"""
struct SpinExpectation{D} <: Observable
    direction::D
    truncate::Bool
end

SpinExpectation(; truncate::Bool = false) = SpinExpectation(nothing, truncate)

function SpinExpectation(direction::AbstractVector; truncate::Bool = false)
    length(direction) == 3 ||
        throw(ArgumentError("a spin-projection direction must have length three"))
    direction_norm = norm(direction)
    isfinite(direction_norm) && !iszero(direction_norm) ||
        throw(ArgumentError("a spin-projection direction must be finite and nonzero"))
    return SpinExpectation(Vec3(direction / direction_norm), truncate)
end

result_name(::SpinExpectation) = :spin_expectation
requirements(::SpinExpectation) = _ObservableRequirements(
    (:hamiltonian, :spin), (:eigensystem, :operator_rotation)
)

function _allocate_observable_result(
        observable::SpinExpectation,
        model::InterpolationModel,
        number_kpoints::Integer,
    )
    T = typeof(real(zero(eltype(model.operators.spin.coefficients))))
    isnothing(observable.direction) &&
        return Array{T}(undef, 3, n_wannier(model), number_kpoints)
    return Matrix{T}(undef, n_wannier(model), number_kpoints)
end

function _spin_expectation!(
        destination::AbstractArray{<:Real},
        observable::SpinExpectation,
        spin::AbstractArray{<:Complex, 4},
        eigenvectors::AbstractArray{<:Complex, 3},
        product::AbstractMatrix,
    )
    number_wannier, _, number_components, number_kpoints = size(spin)
    number_components == 3 ||
        throw(DimensionMismatch("the primitive spin operator must have three components"))
    size(eigenvectors) == (number_wannier, number_wannier, number_kpoints) ||
        throw(DimensionMismatch("spin and eigenvector batches differ"))
    size(product) == (number_wannier, number_wannier) ||
        throw(DimensionMismatch("spin-rotation matrix buffer has the wrong shape"))

    if isnothing(observable.direction)
        size(destination) == (3, number_wannier, number_kpoints) ||
            throw(DimensionMismatch("spin-expectation destination has the wrong shape"))
    else
        size(destination) == (number_wannier, number_kpoints) ||
            throw(DimensionMismatch("projected-spin destination has the wrong shape"))
        fill!(destination, 0)
    end

    for kpoint_index in 1:number_kpoints
        vectors = view(eigenvectors, :, :, kpoint_index)
        for component_index in 1:3
            component = view(spin, :, :, component_index, kpoint_index)
            mul!(product, component, vectors)
            for band_index in 1:number_wannier
                expectation = real(
                    dot(view(vectors, :, band_index), view(product, :, band_index))
                )
                if isnothing(observable.direction)
                    destination[component_index, band_index, kpoint_index] = expectation
                else
                    destination[band_index, kpoint_index] +=
                        observable.direction[component_index] * expectation
                end
            end
        end
    end
    observable.truncate &&
        clamp!(destination, -one(eltype(destination)), one(eltype(destination)))
    return destination
end

function _assemble_observable!(
        destination,
        observable::SpinExpectation,
        intermediates,
        workspace,
    )
    return _spin_expectation!(
        destination,
        observable,
        intermediates.primitive_values.spin,
        intermediates.eigenvectors,
        workspace.operator_product,
    )
end
