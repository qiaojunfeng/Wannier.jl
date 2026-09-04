export BlochOperator, InterpolationModel
export WignerSeitz, MinimumDistance
export Even, Odd, Scalar, PolarVector, AxialVector, CartesianTensor
export component_shape

abstract type TimeReversalParity end

"""Time-reversal-even transformation parity."""
struct Even <: TimeReversalParity end

"""Time-reversal-odd transformation parity."""
struct Odd <: TimeReversalParity end

abstract type OperatorLaw end

"""A scalar operator with the specified time-reversal parity."""
struct Scalar{P <: TimeReversalParity} <: OperatorLaw
    time_reversal::P
end

Scalar(; time_reversal::TimeReversalParity = Even()) = Scalar(time_reversal)

"""A polar-vector operator with the specified time-reversal parity."""
struct PolarVector{P <: TimeReversalParity} <: OperatorLaw
    time_reversal::P
end

PolarVector(; time_reversal::TimeReversalParity = Even()) = PolarVector(time_reversal)

"""An axial-vector operator with the specified time-reversal parity."""
struct AxialVector{P <: TimeReversalParity} <: OperatorLaw
    time_reversal::P
end

AxialVector(; time_reversal::TimeReversalParity = Odd()) = AxialVector(time_reversal)

"""
    CartesianTensor(rank; axial = false, time_reversal = Even())

A rank-`rank` Cartesian tensor operator. Each component axis has length three.
Set `axial = true` when the tensor acquires a determinant factor under improper
spatial operations.
"""
struct CartesianTensor{P <: TimeReversalParity} <: OperatorLaw
    rank::Int
    axial::Bool
    time_reversal::P

    function CartesianTensor(
            rank::Integer,
            axial::Bool,
            time_reversal::P,
        ) where {P <: TimeReversalParity}
        rank >= 0 || throw(ArgumentError("a Cartesian tensor rank must be nonnegative"))
        return new{P}(rank, axial, time_reversal)
    end
end

function CartesianTensor(
        rank::Integer;
        axial::Bool = false,
        time_reversal::TimeReversalParity = Even(),
    )
    return CartesianTensor(rank, axial, time_reversal)
end

component_shape(::Scalar) = ()
component_shape(::Union{PolarVector, AxialVector}) = (3,)
component_shape(law::CartesianTensor) = ntuple(_ -> 3, law.rank)

"""
    BlochOperator(values; law, hermitian)

Matrix elements of a lattice-periodic operator on the source k-point mesh.
Full matrix elements use the layout
`n_bands × n_bands × component_shape... × n_kpoints`. A scalar operator that is
diagonal in the input Bloch basis may instead use `n_bands × n_kpoints`.
"""
struct BlochOperator{A <: AbstractArray, L <: OperatorLaw}
    values::A
    law::L
    hermitian::Bool

    function BlochOperator(
            values::A,
            law::L,
            hermitian::Bool,
        ) where {A <: AbstractArray, L <: OperatorLaw}
        eltype(values) <: Number ||
            throw(ArgumentError("Bloch-operator values must be numeric"))
        ndims(values) >= 2 ||
            throw(ArgumentError("Bloch-operator values must have at least two dimensions"))
        return new{A, L}(values, law, hermitian)
    end
end

function BlochOperator(
        values::AbstractArray{<:Number};
        law::OperatorLaw,
        hermitian::Bool,
    )
    return BlochOperator(values, law, hermitian)
end

abstract type RealSpaceScheme end

"""Select real-space representatives in the Wigner--Seitz cell."""
struct WignerSeitz{T <: Real} <: RealSpaceScheme
    atol::T
    max_cell::Int

    function WignerSeitz(atol::T, max_cell::Integer) where {T <: Real}
        max_cell >= 0 || throw(ArgumentError("max_cell must be nonnegative"))
        return new{T}(atol, Int(max_cell))
    end
end

function WignerSeitz(
        ;
        atol::Real = default_w90_ws_distance_tol(),
        max_cell::Integer = default_w90_ws_search_size(),
    )
    return WignerSeitz(float(atol), Int(max_cell))
end

"""Select the nearest real-space representatives for each Wannier-orbital pair."""
struct MinimumDistance{T <: Real} <: RealSpaceScheme
    atol::T
    max_cell::Int

    function MinimumDistance(atol::T, max_cell::Integer) where {T <: Real}
        max_cell >= 0 || throw(ArgumentError("max_cell must be nonnegative"))
        return new{T}(atol, Int(max_cell))
    end
end

function MinimumDistance(
        ;
        atol::Real = default_w90_ws_distance_tol(),
        max_cell::Integer = default_w90_ws_search_size(),
    )
    return MinimumDistance(float(atol), Int(max_cell))
end

"""
    RealSpaceDomain(lattice, vectors)

A canonical finite set of lattice translations shared by every primitive
operator in an [`InterpolationModel`](@ref). `vectors` are integer fractional
coordinates; `cartesian_vectors` are the corresponding vectors in Cartesian
coordinates.
"""
struct RealSpaceDomain{T <: Real}
    vectors::Vector{Vec3{Int}}
    cartesian_vectors::Vector{Vec3{T}}
    vector_index::Dict{Vec3{Int}, Int}

    function RealSpaceDomain{T}(
            vectors::Vector{Vec3{Int}},
            cartesian_vectors::Vector{Vec3{T}},
            vector_index::Dict{Vec3{Int}, Int},
            ::Val{:validated},
        ) where {T <: Real}
        return new{T}(vectors, cartesian_vectors, vector_index)
    end
end

function RealSpaceDomain(lattice::AbstractMatrix, vectors::AbstractVector)
    size(lattice) == (3, 3) || throw(ArgumentError("lattice must be a 3 × 3 matrix"))
    isempty(vectors) && throw(ArgumentError("a real-space domain cannot be empty"))

    integer_vectors = map(vectors) do vector
        length(vector) == 3 || throw(ArgumentError("each real-space vector must have length three"))
        return Vec3{Int}(vector)
    end
    length(unique(integer_vectors)) == length(integer_vectors) ||
        throw(ArgumentError("real-space vectors must be unique"))
    sort!(integer_vectors; by = Tuple)

    T = float(eltype(lattice))
    lattice_matrix = Mat3{T}(lattice)
    cartesian_vectors = map(integer_vectors) do vector
        Vec3{T}(lattice_matrix * vector)
    end
    vector_index = Dict(vector => index for (index, vector) in enumerate(integer_vectors))
    return RealSpaceDomain{T}(
        integer_vectors, cartesian_vectors, vector_index, Val(:validated)
    )
end

Base.length(real_space::RealSpaceDomain) = length(real_space.vectors)
Base.getindex(real_space::RealSpaceDomain, index::Integer) = real_space.vectors[index]
Base.iterate(real_space::RealSpaceDomain, state...) = iterate(real_space.vectors, state...)
n_Rvectors(real_space::RealSpaceDomain) = length(real_space)

"""
    RealSpaceOperator(coefficients, law, real_space; hermitian)

A primitive Wannier-basis operator on a [`RealSpaceDomain`](@ref). Coefficients
use the layout `n_wannier × n_wannier × component_shape... × n_Rvectors`.
"""
struct RealSpaceOperator{A <: AbstractArray, L <: OperatorLaw}
    coefficients::A
    law::L
    hermitian::Bool

    function RealSpaceOperator(
            coefficients::A,
            law::L,
            hermitian::Bool,
            ::Val{:validated},
        ) where {A <: AbstractArray, L <: OperatorLaw}
        return new{A, L}(coefficients, law, hermitian)
    end
end

function RealSpaceOperator(
        coefficients::AbstractArray{<:Number},
        law::OperatorLaw,
        real_space::RealSpaceDomain,
        ;
        hermitian::Bool,
    )
    ndims(coefficients) >= 3 ||
        throw(ArgumentError("real-space coefficients must have at least three dimensions"))
    size(coefficients, 1) == size(coefficients, 2) ||
        throw(ArgumentError("the two Wannier matrix dimensions must be equal"))
    size(coefficients, ndims(coefficients)) == length(real_space) ||
        throw(ArgumentError("the final coefficient axis must match the real-space domain"))

    actual_shape = Tuple(size(coefficients)[3:(end - 1)])
    expected_shape = component_shape(law)
    actual_shape == expected_shape || throw(
        ArgumentError(
            "operator component shape $actual_shape does not match the law's shape $expected_shape",
        ),
    )
    stored_coefficients = if coefficients isa Array && eltype(coefficients) <: Complex
        coefficients
    else
        complex.(coefficients)
    end
    return RealSpaceOperator(stored_coefficients, law, hermitian, Val(:validated))
end

component_shape(operator::RealSpaceOperator) = Tuple(size(operator.coefficients)[3:(end - 1)])
n_wannier(operator::RealSpaceOperator) = size(operator.coefficients, 1)
n_Rvectors(operator::RealSpaceOperator) = size(operator.coefficients, ndims(operator.coefficients))

struct InterpolationCrystal{T <: Real}
    lattice::Mat3{T}
    atom_positions::Vector{Vec3{T}}
    atom_labels::Vector{String}
end

struct WannierBasis{T <: Real}
    fractional_centers::Vector{Vec3{T}}
end

"""
    InterpolationModel(model; operators = (;), real_space = MinimumDistance())

A persistent set of primitive real-space operators ready for Wannier
interpolation. All operators share `real_space`, so one Fourier phase block can
be reused throughout a calculation.
"""
struct InterpolationModel{C, B, D, O <: NamedTuple, S}
    crystal::C
    basis::B
    real_space::D
    operators::O
    symmetry::S

    function InterpolationModel(
            crystal::C,
            basis::B,
            real_space::D,
            operators::O,
            symmetry::S,
            ::Val{:validated},
        ) where {C, B, D, O <: NamedTuple, S}
        return new{C, B, D, O, S}(crystal, basis, real_space, operators, symmetry)
    end
end

function InterpolationModel(
        crystal,
        basis,
        real_space::RealSpaceDomain,
        operators::NamedTuple,
        symmetry,
    )
    isempty(operators) && throw(ArgumentError("an interpolation model needs an operator"))
    all(operator -> operator isa RealSpaceOperator, values(operators)) ||
        throw(ArgumentError("all interpolation-model operators must be RealSpaceOperator values"))

    first_operator = first(values(operators))
    nwann = n_wannier(first_operator)
    for (name, operator) in pairs(operators)
        n_wannier(operator) == nwann ||
            throw(ArgumentError("operator :$name has a different Wannier dimension"))
        n_Rvectors(operator) == length(real_space) ||
            throw(ArgumentError("operator :$name does not use the common real-space domain"))
    end
    return InterpolationModel(
        crystal, basis, real_space, operators, symmetry, Val(:validated)
    )
end

n_wannier(model::InterpolationModel) = n_wannier(first(values(model.operators)))
n_Rvectors(model::InterpolationModel) = length(model.real_space)
CrystalBase.real_lattice(model::InterpolationModel) = model.crystal.lattice
CrystalBase.reciprocal_lattice(model::InterpolationModel) = reciprocal_lattice(real_lattice(model))
