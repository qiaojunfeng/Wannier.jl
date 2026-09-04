export HamiltonianPosition, PositionHamiltonianPosition

struct _HamiltonianPositionLaw <: OperatorLaw end
struct _PositionHamiltonianPositionLaw <: OperatorLaw end

component_shape(::_HamiltonianPositionLaw) = (3,)
component_shape(::_PositionHamiltonianPositionLaw) = (3, 3)

"""
    HamiltonianPosition()

Construct the Wannier-gauge matrix elements
``\\langle\\mathbf{0}m|H(\\mathbf{r}-\\mathbf{R})|\\mathbf{R}n\\rangle``
needed by [`OrbitalMagnetization`](@ref). Supply the recipe as the
`:hamiltonian_position` entry of `InterpolationModel(...; operators = ...)`.
"""
struct HamiltonianPosition end

_operator_description(::HamiltonianPosition) = (;
    law = _HamiltonianPositionLaw(), hermitian = false,
)

function _validate_operator_input(
        ::Symbol,
        ::HamiltonianPosition,
        ::Integer,
        ::Integer,
    )
    return nothing
end

function _transform_operator_input(::HamiltonianPosition, model::Model)
    kstencil = model.kstencil
    hamiltonian_overlap = compute_hamiltonian_times_position_kspace(
        model.overlaps,
        kstencil.kpb_k,
        kstencil.kpb_G,
        model.eigenvalues,
    )
    wannier_hamiltonian_overlap = transform_gauge(
        hamiltonian_overlap, kstencil.kpb_k, model.gauges
    )

    number_wannier = n_wannier(model)
    number_kpoints = n_kpoints(model)
    T = eltype(wannier_hamiltonian_overlap)
    values = zeros(T, number_wannier, number_wannier, 3, number_kpoints)
    reciprocal = reciprocal_lattice(model)
    for kpoint_index in 1:number_kpoints
        kpoint = kstencil.kpoints[kpoint_index]
        for neighbor_index in 1:n_bvectors(model)
            neighbor_kpoint = kstencil.kpb_k[neighbor_index, kpoint_index]
            reciprocal_shift = kstencil.kpb_G[neighbor_index, kpoint_index]
            neighbor_vector = reciprocal * (
                kstencil.kpoints[neighbor_kpoint] + reciprocal_shift - kpoint
            )
            weight = kstencil.bweights[neighbor_index]
            matrix = view(
                wannier_hamiltonian_overlap,
                :,
                :,
                neighbor_index,
                kpoint_index,
            )
            for component_index in 1:3
                view(values, :, :, component_index, kpoint_index) .+=
                    (im * weight * neighbor_vector[component_index]) .* matrix
            end
        end
    end
    return values
end

"""
    PositionHamiltonianPosition(uHu; force_hermiticity = true)

Construct the Wannier-gauge matrix elements
``\\langle\\mathbf{0}m|r_\\alpha H(r_\\beta-R_\\beta)|\\mathbf{R}n\\rangle``
needed by [`OrbitalMagnetization`](@ref). `uHu` is the neighbor-pair matrix
returned as `read_uHu(path).uHu`. Supply the recipe as the
`:position_hamiltonian_position` entry of `InterpolationModel(...;
operators = ...)`.
"""
struct PositionHamiltonianPosition{U}
    uHu::U
    force_hermiticity::Bool
end

function PositionHamiltonianPosition(
        uHu;
        force_hermiticity::Bool = default_w90_berry_duHdu_force_hermiticity(),
    )
    return PositionHamiltonianPosition(uHu, force_hermiticity)
end

_operator_description(::PositionHamiltonianPosition) = (;
    law = _PositionHamiltonianPositionLaw(), hermitian = false,
)

function _validate_operator_input(
        name::Symbol,
        operator::PositionHamiltonianPosition,
        number_bands::Integer,
        number_kpoints::Integer,
    )
    uHu = operator.uHu
    length(uHu) == number_kpoints || throw(
        DimensionMismatch(
            "operator :$name needs $number_kpoints uHu k points, got $(length(uHu))",
        ),
    )
    for neighbor_pairs in uHu
        size(neighbor_pairs, 1) == size(neighbor_pairs, 2) || throw(
            DimensionMismatch("operator :$name needs square neighbor-pair tables"),
        )
        for matrix in neighbor_pairs
            size(matrix) == (number_bands, number_bands) || throw(
                DimensionMismatch(
                    "operator :$name uHu matrices must have size " *
                        "($number_bands, $number_bands)",
                ),
            )
        end
    end
    return nothing
end

function _hermitize_position_hamiltonian_position!(values::AbstractArray)
    number_wannier = size(values, 1)
    number_kpoints = size(values, 5)
    for kpoint_index in 1:number_kpoints
        for column in 1:number_wannier, row in 1:number_wannier
            for first_component in 1:3, second_component in 1:first_component
                values[row, column, first_component, second_component, kpoint_index] =
                    conj(
                    values[
                        column,
                        row,
                        second_component,
                        first_component,
                        kpoint_index,
                    ],
                )
            end
        end
    end
    return values
end

function _transform_operator_input(
        operator::PositionHamiltonianPosition,
        model::Model,
    )
    kstencil = model.kstencil
    number_neighbors = n_bvectors(model)
    all(size(neighbor_pairs) == (number_neighbors, number_neighbors)
        for neighbor_pairs in operator.uHu) || throw(
        DimensionMismatch(
            "uHu neighbor-pair tables must have size " *
                "($number_neighbors, $number_neighbors)",
        ),
    )

    number_wannier = n_wannier(model)
    number_kpoints = n_kpoints(model)
    sample = first(first(operator.uHu))
    T = promote_type(eltype(sample), eltype(model.gauges))
    values = zeros(T, number_wannier, number_wannier, 3, 3, number_kpoints)
    reciprocal = reciprocal_lattice(model)
    for kpoint_index in 1:number_kpoints
        kpoint = kstencil.kpoints[kpoint_index]
        for second_neighbor in 1:number_neighbors, first_neighbor in 1:number_neighbors
            first_kpoint = kstencil.kpb_k[first_neighbor, kpoint_index]
            second_kpoint = kstencil.kpb_k[second_neighbor, kpoint_index]
            first_shift = kstencil.kpb_G[first_neighbor, kpoint_index]
            second_shift = kstencil.kpb_G[second_neighbor, kpoint_index]
            first_vector = reciprocal * (
                kstencil.kpoints[first_kpoint] + first_shift - kpoint
            )
            second_vector = reciprocal * (
                kstencil.kpoints[second_kpoint] + second_shift - kpoint
            )
            first_gauge = view(model.gauges, :, :, first_kpoint)
            second_gauge = view(model.gauges, :, :, second_kpoint)
            rotated =
                first_gauge' *
                operator.uHu[kpoint_index][first_neighbor, second_neighbor] *
                second_gauge
            weight =
                kstencil.bweights[first_neighbor] *
                kstencil.bweights[second_neighbor]
            for second_component in 1:3, first_component in 1:3
                factor =
                    weight *
                    first_vector[first_component] *
                    second_vector[second_component]
                view(
                    values,
                    :,
                    :,
                    first_component,
                    second_component,
                    kpoint_index,
                ) .+= factor .* rotated
            end
        end
    end
    operator.force_hermiticity &&
        _hermitize_position_hamiltonian_position!(values)
    return values
end
