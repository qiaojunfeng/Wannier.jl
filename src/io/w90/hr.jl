export read_w90_hr, write_w90_hr

"""
    read_w90_hr(prefix, lattice; fractional_centers, atom_positions=[],
                atom_labels=[], symmetry=nothing)

Read `prefix_hr.dat` and its optional `prefix_wsvec.dat` directly into a
Hamiltonian [`InterpolationModel`](@ref). `fractional_centers` are required
because an `hr.dat` file carries no Wannier-center information.
"""
function read_w90_hr(
        prefix::AbstractString,
        lattice::AbstractMatrix;
        fractional_centers::AbstractVector,
        atom_positions::AbstractVector = Vec3{Float64}[],
        atom_labels::AbstractVector = String[],
        symmetry = nothing,
    )
    hrdat = read_w90_hr_dat(prefix * "_hr.dat")
    return InterpolationModel(
        hrdat,
        lattice;
        fractional_centers,
        wsvec = _read_optional_wsvec(prefix),
        atom_positions,
        atom_labels,
        symmetry,
    )
end

"""
    write_w90_hr(prefix, model; write_wsvec=true)

Write the Hamiltonian primitive of an [`InterpolationModel`](@ref) to
`prefix_hr.dat`. The packed real-space domain is written with unit
degeneracies. By default a matching non-MDRS `prefix_wsvec.dat` records that
already-expanded domain.
"""
function write_w90_hr(
        prefix::AbstractString,
        model::InterpolationModel;
        write_wsvec::Bool = true,
    )
    haskey(model.operators, :hamiltonian) ||
        throw(ArgumentError("the interpolation model has no :hamiltonian operator"))
    hamiltonian = model.operators.hamiltonian
    component_shape(hamiltonian) == () ||
        throw(ArgumentError("the :hamiltonian operator must be scalar"))

    T = eltype(real_lattice(model))
    coefficients = Array{Complex{T}}(hamiltonian.coefficients)
    vectors = model.real_space.vectors
    hrdat = WannierIO.HrDat(
        "Written by Wannier.jl",
        vectors,
        ones(Int, length(vectors)),
        coefficients,
    )
    WannierIO.write_w90_hr_dat(prefix * "_hr.dat", hrdat)
    if write_wsvec
        wsvec = WannierIO.WsvecDat(
            "Written by Wannier.jl", vectors, n_wannier(model)
        )
        WannierIO.write_w90_wsvec_dat(prefix * "_wsvec.dat", wsvec)
    end
    return nothing
end
