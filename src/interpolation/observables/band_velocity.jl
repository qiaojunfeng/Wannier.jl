export BandVelocity

"""
    BandVelocity()

Request the Cartesian band-energy derivative
``\\hbar \\mathbf{v}_n = \\partial \\varepsilon_n / \\partial \\mathbf{k}``.
The result is stored in eV Å with layout
`3 × n_wannier × n_kpoints`. At an exact degeneracy, individual diagonal
velocities depend on the eigenvector basis within the degenerate subspace; a
full velocity matrix is required when that subspace information is needed.
"""
struct BandVelocity <: Observable end

result_name(::BandVelocity) = :band_velocity
requirements(::BandVelocity) = _ObservableRequirements(
    (:hamiltonian,), (:eigensystem, :hamiltonian_gradient)
)

function _allocate_observable_result(
        ::BandVelocity,
        model::InterpolationModel,
        number_kpoints::Integer,
    )
    T = typeof(real(zero(eltype(model.operators.hamiltonian.coefficients))))
    return Array{T}(undef, 3, n_wannier(model), number_kpoints)
end

function _band_velocity!(
        destination::AbstractArray{<:Real, 3},
        hamiltonian_gradient::AbstractArray{<:Complex, 4},
        eigenvectors::AbstractArray{<:Complex, 3},
        product::AbstractMatrix,
    )
    number_wannier, _, _, number_kpoints = size(hamiltonian_gradient)
    size(destination) == (3, number_wannier, number_kpoints) ||
        throw(DimensionMismatch("band-velocity destination has the wrong shape"))
    size(eigenvectors) == (number_wannier, number_wannier, number_kpoints) ||
        throw(DimensionMismatch("Hamiltonian-gradient and eigenvector batches differ"))
    size(product) == (number_wannier, number_wannier) ||
        throw(DimensionMismatch("band-velocity matrix buffer has the wrong shape"))

    for kpoint_index in 1:number_kpoints
        vectors = view(eigenvectors, :, :, kpoint_index)
        for cartesian_component in 1:3
            derivative = view(
                hamiltonian_gradient, :, :, cartesian_component, kpoint_index
            )
            mul!(product, derivative, vectors)
            for band_index in 1:number_wannier
                destination[cartesian_component, band_index, kpoint_index] = real(
                    dot(view(vectors, :, band_index), view(product, :, band_index))
                )
            end
        end
    end
    return destination
end

function _assemble_observable!(
        destination::AbstractArray{<:Real, 3},
        ::BandVelocity,
        intermediates,
        workspace,
    )
    return _band_velocity!(
        destination,
        intermediates.hamiltonian_gradient,
        intermediates.eigenvectors,
        workspace.operator_product,
    )
end
