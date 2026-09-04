export read_w90_tb, write_w90_tb

function _read_optional_wsvec(prefix::AbstractString)
    path = prefix * "_wsvec.dat"
    return isfile(path) ? read_w90_wsvec_dat(path) : nothing
end

"""
    read_w90_tb(prefix; fractional_centers=nothing, atom_positions=[],
                atom_labels=[], symmetry=nothing)

Read `prefix_tb.dat` and its optional `prefix_wsvec.dat` directly into an
[`InterpolationModel`](@ref). The model contains both `:hamiltonian` and
`:berry_connection` primitive operators. If `fractional_centers` is omitted,
the centers are inferred from the diagonal `R = 0` position matrix elements.
"""
function read_w90_tb(
        prefix::AbstractString;
        fractional_centers::Union{Nothing, AbstractVector} = nothing,
        atom_positions::AbstractVector = Vec3{Float64}[],
        atom_labels::AbstractVector = String[],
        symmetry = nothing,
    )
    tbdat = read_w90_tb_dat(prefix * "_tb.dat")
    return InterpolationModel(
        tbdat;
        fractional_centers,
        wsvec = _read_optional_wsvec(prefix),
        atom_positions,
        atom_labels,
        symmetry,
    )
end

"""
    write_w90_tb(prefix, model; write_wsvec=true)

Write the Hamiltonian and Berry-connection primitives of an
[`InterpolationModel`](@ref) to `prefix_tb.dat`. The packed real-space domain
is written with unit degeneracies. By default a matching non-MDRS
`prefix_wsvec.dat` records that already-expanded domain.
"""
function write_w90_tb(
        prefix::AbstractString,
        model::InterpolationModel;
        write_wsvec::Bool = true,
    )
    haskey(model.operators, :hamiltonian) ||
        throw(ArgumentError("the interpolation model has no :hamiltonian operator"))
    haskey(model.operators, :berry_connection) || throw(
        ArgumentError("the interpolation model has no :berry_connection operator"),
    )
    hamiltonian = model.operators.hamiltonian
    position = model.operators.berry_connection
    component_shape(hamiltonian) == () ||
        throw(ArgumentError("the :hamiltonian operator must be scalar"))
    component_shape(position) == (3,) ||
        throw(ArgumentError("the :berry_connection operator must have three components"))

    T = eltype(real_lattice(model))
    complex_type = Complex{T}
    hamiltonian_coefficients = Array{complex_type}(hamiltonian.coefficients)
    position_coefficients = Array{complex_type}(position.coefficients)
    vectors = model.real_space.vectors
    degeneracies = ones(Int, length(vectors))
    tbdat = WannierIO.TbDat(
        "Written by Wannier.jl",
        real_lattice(model),
        vectors,
        degeneracies,
        hamiltonian_coefficients,
        Array(view(position_coefficients, :, :, 1, :)),
        Array(view(position_coefficients, :, :, 2, :)),
        Array(view(position_coefficients, :, :, 3, :)),
    )
    WannierIO.write_w90_tb_dat(prefix * "_tb.dat", tbdat)
    if write_wsvec
        wsvec = WannierIO.WsvecDat(
            "Written by Wannier.jl", vectors, n_wannier(model)
        )
        WannierIO.write_w90_wsvec_dat(prefix * "_wsvec.dat", wsvec)
    end
    return nothing
end
