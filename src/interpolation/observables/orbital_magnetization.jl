export OrbitalMagnetization

"""
    OrbitalMagnetization(fermi_energy; degen_tol = 1e-4)

Request the antisymmetric Cartesian orbital-magnetization integrand of LVTS12,
with shape `3 × 3 × n_kpoints` and units eV Å². The model must contain
`:berry_connection`, `:hamiltonian_position`, and
`:position_hamiltonian_position` primitives constructed with
[`BerryConnection`](@ref), [`HamiltonianPosition`](@ref), and
[`PositionHamiltonianPosition`](@ref), respectively.

`degen_tol` is reserved for the degenerate-perturbation formulation; the
current recipe follows Wannier90's default non-perturbative branch.
"""
struct OrbitalMagnetization{E} <: Observable
    fermi_energy::E
    degen_tol::Float64
end

function OrbitalMagnetization(
        fermi_energy::Real;
        degen_tol::Real = default_w90_berry_degen_tol(),
    )
    return OrbitalMagnetization(float(fermi_energy), float(degen_tol))
end

result_name(::OrbitalMagnetization) = :orbital_magnetization
requirements(::OrbitalMagnetization) = _ObservableRequirements(
    (
        :hamiltonian,
        :berry_connection,
        :hamiltonian_position,
        :position_hamiltonian_position,
    ),
    (
        :eigensystem,
        :hamiltonian_gradient,
        :berry_connection_gradient,
        :operator_rotation,
    ),
)

function _allocate_observable_result(
        ::OrbitalMagnetization,
        model::InterpolationModel,
        number_kpoints::Integer,
    )
    T = typeof(real(zero(eltype(model.operators.hamiltonian.coefficients))))
    return Array{T}(undef, 3, 3, number_kpoints)
end

struct _OrbitalMagnetizationWorkspace{B, T <: Complex}
    berry::B
    hamiltonian_position::Array{T, 3}
    lambda::Array{T, 4}
end

function _allocate_observable_workspace(
        ::OrbitalMagnetization,
        model::InterpolationModel,
    )
    number_wannier = n_wannier(model)
    berry = _allocate_berry_curvature_workspace(model)
    T = promote_type(
        eltype(model.operators.hamiltonian_position.coefficients),
        eltype(model.operators.position_hamiltonian_position.coefficients),
        eltype(berry.connection),
    )
    return _OrbitalMagnetizationWorkspace(
        berry,
        Array{T}(undef, number_wannier, number_wannier, 3),
        Array{T}(undef, number_wannier, number_wannier, 3, 3),
    )
end

function _prepare_orbital_magnetization_workspace!(
        workspace::_OrbitalMagnetizationWorkspace,
        intermediates,
        kpoint_index::Integer,
    )
    berry = workspace.berry
    _prepare_berry_workspace!(berry, intermediates, kpoint_index)
    eigenvectors = view(intermediates.eigenvectors, :, :, kpoint_index)
    hamiltonian_position = intermediates.primitive_values.hamiltonian_position
    position_hamiltonian_position =
        intermediates.primitive_values.position_hamiltonian_position

    for component_index in 1:3
        _rotate_operator_component!(
            view(workspace.hamiltonian_position, :, :, component_index),
            view(
                hamiltonian_position,
                :,
                :,
                component_index,
                kpoint_index,
            ),
            eigenvectors,
            berry.product,
        )
    end
    for second_component in 1:3, first_component in 1:3
        for column in axes(berry.source, 2), row in axes(berry.source, 1)
            berry.source[row, column] = im * (
                position_hamiltonian_position[
                    row,
                    column,
                    first_component,
                    second_component,
                    kpoint_index,
                ] -
                    conj(
                    position_hamiltonian_position[
                        column,
                        row,
                        first_component,
                        second_component,
                        kpoint_index,
                    ],
                )
            )
        end
        _rotate_operator_component!(
            view(
                workspace.lambda,
                :,
                :,
                first_component,
                second_component,
            ),
            berry.source,
            eigenvectors,
            berry.product,
        )
    end
    return workspace
end

function _orbital_magnetization!(
        destination::AbstractArray{<:Real, 3},
        observable::OrbitalMagnetization,
        intermediates,
        workspace::_OrbitalMagnetizationWorkspace,
    )
    number_wannier, number_kpoints = size(intermediates.eigenvalues)
    size(destination) == (3, 3, number_kpoints) || throw(
        DimensionMismatch("orbital-magnetization destination has the wrong shape"),
    )

    for kpoint_index in 1:number_kpoints
        _prepare_orbital_magnetization_workspace!(
            workspace, intermediates, kpoint_index
        )
        berry = workspace.berry
        eigenvalues = view(intermediates.eigenvalues, :, kpoint_index)
        for band_index in eachindex(eigenvalues)
            berry.occupation[band_index] =
                eigenvalues[band_index] <= observable.fermi_energy
        end

        connection = berry.connection
        curvature = berry.curvature
        gauge_derivative = berry.gauge_derivative
        hamiltonian_position = workspace.hamiltonian_position
        lambda = workspace.lambda
        for second_component in 1:3, first_component in 1:3
            curvature_trace = zero(eltype(connection))
            energy_curvature_trace = zero(eltype(connection))
            lambda_trace = zero(eltype(connection))
            connection_pair = zero(eltype(connection))
            mixed_pair = zero(eltype(connection))
            energy_mixed_pair = zero(eltype(connection))
            hamiltonian_connection_pair = zero(eltype(connection))

            for occupied_band in 1:number_wannier
                berry.occupation[occupied_band] || continue
                occupied_energy = eigenvalues[occupied_band]
                curvature_element = curvature[
                    occupied_band,
                    occupied_band,
                    first_component,
                    second_component,
                ]
                curvature_trace += curvature_element
                energy_curvature_trace += occupied_energy * curvature_element
                lambda_trace += lambda[
                    occupied_band,
                    occupied_band,
                    first_component,
                    second_component,
                ]

                for inner_band in 1:number_wannier
                    if berry.occupation[inner_band]
                        connection_pair +=
                            occupied_energy *
                            connection[
                                occupied_band,
                                inner_band,
                                first_component,
                            ] *
                            connection[
                                inner_band,
                                occupied_band,
                                second_component,
                            ]
                        continue
                    end

                    j_minus_first =
                        im *
                        gauge_derivative[
                            occupied_band,
                            inner_band,
                            first_component,
                        ]
                    j_minus_second =
                        im *
                        gauge_derivative[
                            occupied_band,
                            inner_band,
                            second_component,
                        ]
                    j_plus_second =
                        im *
                        gauge_derivative[
                            inner_band,
                            occupied_band,
                            second_component,
                        ]
                    mixed_element =
                        connection[
                        occupied_band,
                        inner_band,
                        first_component,
                    ] * j_plus_second +
                        j_minus_first * (
                        connection[
                            inner_band,
                            occupied_band,
                            second_component,
                        ] + j_plus_second
                    )
                    mixed_pair += mixed_element
                    energy_mixed_pair += occupied_energy * mixed_element
                    hamiltonian_connection_pair +=
                        j_minus_first *
                        hamiltonian_position[
                            inner_band,
                            occupied_band,
                            second_component,
                        ] -
                        j_minus_second *
                        hamiltonian_position[
                            inner_band,
                            occupied_band,
                            first_component,
                        ] +
                        j_minus_first *
                        eigenvalues[inner_band] *
                        j_plus_second
                end
            end

            berry_curvature = real(curvature_trace) - 2imag(mixed_pair)
            itinerant =
                real(energy_curvature_trace) + 2imag(connection_pair) -
                2imag(energy_mixed_pair)
            local_circulation =
                real(lambda_trace) - 2imag(connection_pair) -
                2imag(hamiltonian_connection_pair)
            destination[first_component, second_component, kpoint_index] =
                itinerant + local_circulation -
                2observable.fermi_energy * berry_curvature
        end
    end
    return destination
end

function _assemble_observable!(
        destination,
        observable::OrbitalMagnetization,
        intermediates,
        workspace,
    )
    return _orbital_magnetization!(
        destination,
        observable,
        intermediates,
        workspace.observable_workspaces.orbital_magnetization,
    )
end
