export BerryConnection

struct _BerryConnectionLaw <: OperatorLaw end

component_shape(::_BerryConnectionLaw) = (3,)

"""
    BerryConnection(; imlog_diag = true, force_hermiticity = true)

Construction recipe for the finite-difference Wannier-gauge Berry connection.
Supply it as the `:berry_connection` entry of `InterpolationModel(...;
operators = ...)`. `imlog_diag` and `force_hermiticity` have the same meaning as
in [`compute_berry_connection_kspace`](@ref).

The connection has an inhomogeneous gauge-transformation law and is therefore
kept distinct from a pointwise [`BlochOperator`](@ref).
"""
struct BerryConnection
    imlog_diag::Bool
    force_hermiticity::Bool
end

function BerryConnection(
        ;
        imlog_diag::Bool = true,
        force_hermiticity::Bool = default_w90_berry_position_force_hermiticity(),
    )
    return BerryConnection(imlog_diag, force_hermiticity)
end

_operator_description(::BerryConnection) = (;
    law = _BerryConnectionLaw(), hermitian = true,
)

function _validate_operator_input(
        ::Symbol,
        ::BerryConnection,
        ::Integer,
        ::Integer,
    )
    return nothing
end

function _dense_berry_connection(connection::AbstractVector)
    isempty(connection) && throw(ArgumentError("Berry connection cannot be empty"))
    number_kpoints = length(connection)
    number_wannier = size(first(connection), 1)
    size(first(connection), 2) == number_wannier ||
        throw(DimensionMismatch("Berry-connection matrices must be square"))
    T = eltype(eltype(first(connection)))
    values = Array{T}(undef, number_wannier, number_wannier, 3, number_kpoints)
    for kpoint_index in eachindex(connection)
        matrix = connection[kpoint_index]
        size(matrix) == (number_wannier, number_wannier) ||
            throw(DimensionMismatch("Berry-connection matrix sizes differ"))
        for column in 1:number_wannier, row in 1:number_wannier
            value = matrix[row, column]
            length(value) == 3 ||
                throw(DimensionMismatch("Berry connections need three components"))
            for component_index in 1:3
                values[row, column, component_index, kpoint_index] =
                    value[component_index]
            end
        end
    end
    return values
end

function _transform_operator_input(connection::BerryConnection, model::Model)
    values = compute_berry_connection_kspace(
        model,
        model.gauges;
        imlog_diag = connection.imlog_diag,
        force_hermiticity = connection.force_hermiticity,
    )
    return _dense_berry_connection(values)
end
