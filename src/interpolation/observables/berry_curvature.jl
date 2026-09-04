export BerryCurvature, WYSV06, WYSV06BandResolved, LVTS12

"""
    BerryCurvature(fermi_energy; formulation = LVTS12(), degen_tol = 1e-4)
    BerryCurvature(; formulation = WYSV06BandResolved(), degen_tol = 1e-4)

Request the antisymmetric Cartesian Berry-curvature tensor in Å². The
occupied-manifold formulations return `result.berry_curvature` with layout
`3 × 3 × n_kpoints`; `WYSV06BandResolved()` returns
`3 × 3 × n_wannier × n_kpoints`.

The interpolation model must contain a `:berry_connection` primitive, normally
constructed with [`BerryConnection`](@ref). `degen_tol` is reserved for the
degenerate-perturbation formulation; the current recipe follows Wannier90's
default non-perturbative branch.
"""
struct BerryCurvature{F <: AbstractBerryCurvatureInterpolationAlgorithm, E} <: Observable
    formulation::F
    fermi_energy::E
    degen_tol::Float64
end

function BerryCurvature(
        fermi_energy::Real;
        formulation::AbstractBerryCurvatureInterpolationAlgorithm = LVTS12(),
        degen_tol::Real = default_w90_berry_degen_tol(),
    )
    formulation isa WYSV06BandResolved && throw(
        ArgumentError("the band-resolved formulation does not use a Fermi energy"),
    )
    return BerryCurvature(formulation, float(fermi_energy), float(degen_tol))
end

function BerryCurvature(
        ;
        formulation::AbstractBerryCurvatureInterpolationAlgorithm = WYSV06BandResolved(),
        degen_tol::Real = default_w90_berry_degen_tol(),
    )
    formulation isa WYSV06BandResolved || throw(
        ArgumentError("an occupied-manifold formulation requires a Fermi energy"),
    )
    return BerryCurvature(formulation, nothing, float(degen_tol))
end

result_name(::BerryCurvature) = :berry_curvature
requirements(::BerryCurvature) = _ObservableRequirements(
    (:hamiltonian, :berry_connection),
    (
        :eigensystem,
        :hamiltonian_gradient,
        :berry_connection_gradient,
        :operator_rotation,
    ),
)

function _allocate_observable_result(
        observable::BerryCurvature,
        model::InterpolationModel,
        number_kpoints::Integer,
    )
    T = typeof(real(zero(eltype(model.operators.hamiltonian.coefficients))))
    observable.formulation isa WYSV06BandResolved &&
        return Array{T}(undef, 3, 3, n_wannier(model), number_kpoints)
    return Array{T}(undef, 3, 3, number_kpoints)
end

struct _BerryCurvatureWorkspace{T <: Complex}
    hamiltonian_gradient::Array{T, 3}
    connection::Array{T, 3}
    curvature_wannier::Array{T, 4}
    curvature::Array{T, 4}
    gauge_derivative::Array{T, 3}
    occupation::BitVector
    occupation_wannier::Matrix{T}
    gauge_connection::Array{T, 3}
    gauge_connection_minus::Array{T, 3}
    gauge_connection_plus::Array{T, 3}
    source::Matrix{T}
    product::Matrix{T}
end

function _allocate_observable_workspace(
        ::BerryCurvature,
        model::InterpolationModel,
    )
    number_wannier = n_wannier(model)
    T = promote_type(
        eltype(model.operators.hamiltonian.coefficients),
        eltype(model.operators.berry_connection.coefficients),
    )
    return _BerryCurvatureWorkspace(
        Array{T}(undef, number_wannier, number_wannier, 3),
        Array{T}(undef, number_wannier, number_wannier, 3),
        Array{T}(undef, number_wannier, number_wannier, 3, 3),
        Array{T}(undef, number_wannier, number_wannier, 3, 3),
        Array{T}(undef, number_wannier, number_wannier, 3),
        falses(number_wannier),
        Matrix{T}(undef, number_wannier, number_wannier),
        Array{T}(undef, number_wannier, number_wannier, 3),
        Array{T}(undef, number_wannier, number_wannier, 3),
        Array{T}(undef, number_wannier, number_wannier, 3),
        Matrix{T}(undef, number_wannier, number_wannier),
        Matrix{T}(undef, number_wannier, number_wannier),
    )
end

function _rotate_from_band_basis!(
        destination::AbstractMatrix,
        source::AbstractMatrix,
        eigenvectors::AbstractMatrix,
        product::AbstractMatrix,
    )
    mul!(product, source, eigenvectors')
    mul!(destination, eigenvectors, product)
    return destination
end

function _rotate_operator_component!(
        destination::AbstractMatrix,
        source::AbstractMatrix,
        eigenvectors::AbstractMatrix,
        product::AbstractMatrix,
    )
    mul!(product, source, eigenvectors)
    mul!(destination, eigenvectors', product)
    return destination
end

function _prepare_berry_workspace!(
        workspace::_BerryCurvatureWorkspace,
        intermediates,
        kpoint_index::Integer,
    )
    eigenvectors = view(intermediates.eigenvectors, :, :, kpoint_index)
    hamiltonian_gradient = intermediates.primitive_derivatives.hamiltonian
    connection = intermediates.primitive_values.berry_connection
    connection_gradient = intermediates.primitive_derivatives.berry_connection

    for component_index in 1:3
        _rotate_operator_component!(
            view(workspace.hamiltonian_gradient, :, :, component_index),
            view(hamiltonian_gradient, :, :, component_index, kpoint_index),
            eigenvectors,
            workspace.product,
        )
        _rotate_operator_component!(
            view(workspace.connection, :, :, component_index),
            view(connection, :, :, component_index, kpoint_index),
            eigenvectors,
            workspace.product,
        )
    end

    for second_component in 1:3, first_component in 1:3
        curvature_wannier = view(
            workspace.curvature_wannier, :, :, first_component, second_component
        )
        for column in axes(workspace.source, 2), row in axes(workspace.source, 1)
            curvature_wannier[row, column] =
                connection_gradient[
                row, column, second_component, first_component, kpoint_index,
            ] -
                connection_gradient[
                row, column, first_component, second_component, kpoint_index,
            ]
        end
        _rotate_operator_component!(
            view(workspace.curvature, :, :, first_component, second_component),
            curvature_wannier,
            eigenvectors,
            workspace.product,
        )
    end

    eigenvalues = view(intermediates.eigenvalues, :, kpoint_index)
    for component_index in 1:3
        for column in eachindex(eigenvalues), row in eachindex(eigenvalues)
            workspace.gauge_derivative[row, column, component_index] = if row == column
                zero(eltype(workspace.gauge_derivative))
            else
                workspace.hamiltonian_gradient[row, column, component_index] /
                    (eigenvalues[column] - eigenvalues[row])
            end
        end
    end
    return workspace
end


function _prepare_lvts_workspace!(
        workspace::_BerryCurvatureWorkspace,
        intermediates,
        kpoint_index::Integer,
        fermi_energy::Real,
    )
    eigenvalues = view(intermediates.eigenvalues, :, kpoint_index)
    eigenvectors = view(intermediates.eigenvectors, :, :, kpoint_index)
    for band_index in eachindex(eigenvalues)
        workspace.occupation[band_index] = eigenvalues[band_index] <= fermi_energy
    end

    fill!(workspace.source, 0)
    for band_index in eachindex(eigenvalues)
        workspace.source[band_index, band_index] = workspace.occupation[band_index]
    end
    _rotate_from_band_basis!(
        workspace.occupation_wannier,
        workspace.source,
        eigenvectors,
        workspace.product,
    )

    for component_index in 1:3
        for column in eachindex(eigenvalues), row in eachindex(eigenvalues)
            value = im * workspace.gauge_derivative[row, column, component_index]
            workspace.source[row, column] = value
        end
        _rotate_from_band_basis!(
            view(workspace.gauge_connection, :, :, component_index),
            workspace.source,
            eigenvectors,
            workspace.product,
        )

        for column in eachindex(eigenvalues), row in eachindex(eigenvalues)
            value = im * workspace.gauge_derivative[row, column, component_index]
            workspace.source[row, column] =
                workspace.occupation[row] && !workspace.occupation[column] ? value : 0
        end
        _rotate_from_band_basis!(
            view(workspace.gauge_connection_minus, :, :, component_index),
            workspace.source,
            eigenvectors,
            workspace.product,
        )

        for column in eachindex(eigenvalues), row in eachindex(eigenvalues)
            value = im * workspace.gauge_derivative[row, column, component_index]
            workspace.source[row, column] =
                !workspace.occupation[row] && workspace.occupation[column] ? value : 0
        end
        _rotate_from_band_basis!(
            view(workspace.gauge_connection_plus, :, :, component_index),
            workspace.source,
            eigenvectors,
            workspace.product,
        )
    end
    return workspace
end

function _band_resolved_berry_curvature!(
        destination::AbstractArray{<:Real, 4},
        intermediates,
        workspace::_BerryCurvatureWorkspace,
    )
    number_wannier, number_kpoints = size(intermediates.eigenvalues)
    size(destination) == (3, 3, number_wannier, number_kpoints) ||
        throw(DimensionMismatch("band-resolved Berry-curvature destination has the wrong shape"))

    for kpoint_index in 1:number_kpoints
        _prepare_berry_workspace!(workspace, intermediates, kpoint_index)
        A = workspace.connection
        D = workspace.gauge_derivative
        omega = workspace.curvature
        for band_index in 1:number_wannier
            for second_component in 1:3, first_component in 1:3
                DA = zero(eltype(A))
                DA_transpose = zero(eltype(A))
                DD = zero(eltype(A))
                DD_transpose = zero(eltype(A))
                for inner_index in 1:number_wannier
                    DA +=
                        D[band_index, inner_index, first_component] *
                        A[inner_index, band_index, second_component] -
                        A[band_index, inner_index, second_component] *
                        D[inner_index, band_index, first_component]
                    DA_transpose +=
                        D[band_index, inner_index, second_component] *
                        A[inner_index, band_index, first_component] -
                        A[band_index, inner_index, first_component] *
                        D[inner_index, band_index, second_component]
                    DD +=
                        D[band_index, inner_index, first_component] *
                        D[inner_index, band_index, second_component]
                    DD_transpose +=
                        D[band_index, inner_index, second_component] *
                        D[inner_index, band_index, first_component]
                end
                destination[
                    first_component, second_component, band_index, kpoint_index,
                ] = real(
                    omega[
                        band_index, band_index, first_component, second_component,
                    ] -
                        (DA - DA_transpose) - im * (DD - DD_transpose),
                )
            end
        end
    end
    return destination
end

function _occupied_berry_curvature!(
        destination::AbstractArray{<:Real, 3},
        observable::BerryCurvature{<:WYSV06},
        intermediates,
        workspace::_BerryCurvatureWorkspace,
    )
    number_wannier, number_kpoints = size(intermediates.eigenvalues)
    size(destination) == (3, 3, number_kpoints) ||
        throw(DimensionMismatch("Berry-curvature destination has the wrong shape"))

    for kpoint_index in 1:number_kpoints
        _prepare_berry_workspace!(workspace, intermediates, kpoint_index)
        eigenvalues = view(intermediates.eigenvalues, :, kpoint_index)
        A = workspace.connection
        D = workspace.gauge_derivative
        omega = workspace.curvature
        for second_component in 1:3, first_component in 1:3
            value = zero(eltype(A))
            for band_index in 1:number_wannier
                occupation = eigenvalues[band_index] <= observable.fermi_energy
                if occupation
                    value += omega[
                        band_index, band_index, first_component, second_component,
                    ]
                end
                DA = zero(eltype(A))
                DA_transpose = zero(eltype(A))
                DD = zero(eltype(A))
                for inner_index in 1:number_wannier
                    inner_occupation =
                        eigenvalues[inner_index] <= observable.fermi_energy
                    occupation_difference = inner_occupation - occupation
                    DA +=
                        occupation_difference *
                        D[band_index, inner_index, first_component] *
                        A[inner_index, band_index, second_component]
                    DA_transpose +=
                        occupation_difference *
                        D[band_index, inner_index, second_component] *
                        A[inner_index, band_index, first_component]
                    DD +=
                        im * occupation_difference *
                        D[band_index, inner_index, first_component] *
                        D[inner_index, band_index, second_component]
                end
                value += DA - DA_transpose + DD
            end
            destination[first_component, second_component, kpoint_index] = real(value)
        end
    end
    return destination
end

function _occupied_berry_curvature!(
        destination::AbstractArray{<:Real, 3},
        observable::BerryCurvature{<:LVTS12},
        intermediates,
        workspace::_BerryCurvatureWorkspace,
    )
    number_wannier, number_kpoints = size(intermediates.eigenvalues)
    size(destination) == (3, 3, number_kpoints) ||
        throw(DimensionMismatch("Berry-curvature destination has the wrong shape"))

    for kpoint_index in 1:number_kpoints
        _prepare_berry_workspace!(workspace, intermediates, kpoint_index)
        _prepare_lvts_workspace!(
            workspace, intermediates, kpoint_index, observable.fermi_energy
        )
        connection_wannier = view(
            intermediates.primitive_values.berry_connection, :, :, :, kpoint_index
        )
        for second_component in 1:3, first_component in 1:3
            value = zero(eltype(workspace.source))
            correction = zero(eltype(workspace.source))
            for column in 1:number_wannier, row in 1:number_wannier
                value +=
                    workspace.occupation_wannier[row, column] *
                    workspace.curvature_wannier[
                        column, row, first_component, second_component,
                    ]
                correction +=
                    connection_wannier[row, column, first_component] *
                    workspace.gauge_connection_plus[
                        column, row, second_component,
                    ] +
                    workspace.gauge_connection_minus[
                        row, column, first_component,
                    ] *
                    (
                        connection_wannier[column, row, second_component] +
                            workspace.gauge_connection[
                                column, row, second_component,
                            ]
                    )
            end
            destination[first_component, second_component, kpoint_index] =
                real(value) - 2imag(correction)
        end
    end
    return destination
end

function _assemble_observable!(
        destination,
        observable::BerryCurvature,
        intermediates,
        workspace,
    )
    berry_workspace = workspace.observable_workspaces.berry_curvature
    observable.formulation isa WYSV06BandResolved &&
        return _band_resolved_berry_curvature!(
        destination, intermediates, berry_workspace
    )
    observable.formulation isa WYSV06 &&
        return _occupied_berry_curvature!(
        destination, observable, intermediates, berry_workspace
        )
    observable.formulation isa LVTS12 &&
        return _occupied_berry_curvature!(
            destination, observable, intermediates, berry_workspace
        )
    error("unsupported Berry-curvature formulation $(typeof(observable.formulation))")
end
