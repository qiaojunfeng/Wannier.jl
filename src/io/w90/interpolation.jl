function _w90_real_space_selection(
        vectors,
        degeneracies,
        number_wannier::Integer,
        wsvec::Union{Nothing, WannierIO.WsvecDat},
    )
    wigner_seitz = _WignerSeitzSelection(vectors, degeneracies)
    isnothing(wsvec) && return wigner_seitz

    wsvec.Rvectors == wigner_seitz.quotient_vectors || throw(
        ArgumentError("the Hamiltonian and wsvec files use different R vectors"),
    )
    WannierIO.n_wannier(wsvec) == number_wannier || throw(
        DimensionMismatch("the Hamiltonian and wsvec files use different Wannier dimensions"),
    )
    wsvec.mdrs || return wigner_seitz
    isnothing(wsvec.Tvectors) &&
        throw(ArgumentError("an MDRS wsvec file must contain translation vectors"))
    isnothing(wsvec.Tdegens) &&
        throw(ArgumentError("an MDRS wsvec file must contain translation degeneracies"))
    return _minimum_distance_selection(
        wigner_seitz, wsvec.Tvectors, wsvec.Tdegens
    )
end

function _interpolation_model_from_w90_hamiltonian(
        hamiltonian,
        lattice,
        vectors,
        degeneracies;
        fractional_centers,
        wsvec,
        atom_positions,
        atom_labels,
        symmetry,
    )
    isnothing(symmetry) || symmetry isa WannierSymmetry ||
        throw(ArgumentError("symmetry must be a WannierSymmetry or nothing"))
    number_wannier = size(hamiltonian, 1)
    size(hamiltonian, 2) == number_wannier ||
        throw(DimensionMismatch("the Hamiltonian matrix dimensions must be equal"))
    size(hamiltonian, 3) == length(vectors) ||
        throw(DimensionMismatch("Hamiltonian and R-vector counts differ"))

    selection = _w90_real_space_selection(
        vectors, degeneracies, number_wannier, wsvec
    )
    selected_hamiltonian = _apply_real_space_selection(selection, hamiltonian)
    descriptions = (;
        hamiltonian = (;
            law = Scalar(time_reversal = Even()),
            hermitian = true,
        ),
    )
    representative_vectors = _representative_vectors(selection)
    selected_coefficients = (selected_hamiltonian,)
    if !isnothing(symmetry)
        n_wannier(symmetry) == number_wannier || throw(
            DimensionMismatch("Hamiltonian and Wannier symmetry have different basis sizes"),
        )
        representative_vectors, selected_coefficients = _close_real_space_operators(
            lattice, representative_vectors, descriptions, selected_coefficients, symmetry
        )
    end
    domain, operators = _pack_real_space_operators(
        lattice, representative_vectors, descriptions, selected_coefficients
    )
    crystal = _interpolation_crystal(lattice, atom_positions, atom_labels)
    basis = _wannier_basis(
        fractional_centers, number_wannier, eltype(crystal.lattice)
    )
    return InterpolationModel(crystal, basis, domain, operators, symmetry)
end

"""
    InterpolationModel(hrdat, lattice; fractional_centers, wsvec=nothing,
                       atom_positions=[], atom_labels=[])

Construct a Hamiltonian interpolation model directly from Wannier90 `HrDat`
file data. `fractional_centers` are the Wannier centers in lattice coordinates.
When the corresponding `WsvecDat` is supplied, its minimum-distance replicas
are absorbed into the stored coefficients; otherwise the `hrdat` Wigner--Seitz
degeneracies are used.
"""
function InterpolationModel(
        hrdat::WannierIO.HrDat,
        lattice::AbstractMatrix;
        fractional_centers::AbstractVector,
        wsvec::Union{Nothing, WannierIO.WsvecDat} = nothing,
        atom_positions::AbstractVector = Vec3{Float64}[],
        atom_labels::AbstractVector = String[],
        symmetry = nothing,
    )
    return _interpolation_model_from_w90_hamiltonian(
        hrdat.H,
        lattice,
        hrdat.Rvectors,
        hrdat.Rdegens;
        fractional_centers,
        wsvec,
        atom_positions,
        atom_labels,
        symmetry,
    )
end

"""
    BlochOperator(spn::WannierIO.Spn)

Construct the Hermitian axial-vector, time-reversal-odd Bloch operator stored
in a Wannier90 `spn` file. Its logical layout is
`n_bands × n_bands × 3 × n_kpoints`.
"""
function BlochOperator(spn::WannierIO.Spn)
    number_bands, _, number_kpoints = size(spn.Sx)
    values = Array{eltype(spn.Sx)}(
        undef, number_bands, number_bands, 3, number_kpoints
    )
    view(values, :, :, 1, :) .= spn.Sx
    view(values, :, :, 2, :) .= spn.Sy
    view(values, :, :, 3, :) .= spn.Sz
    return BlochOperator(
        values; law = AxialVector(time_reversal = Odd()), hermitian = true
    )
end

"""
    InterpolationModel(tbdat; fractional_centers, wsvec=nothing,
                       atom_positions=[], atom_labels=[])

Construct a Hamiltonian interpolation model directly from Wannier90 `TbDat`
file data. The lattice is read from `tbdat`; other behavior matches the
[`HrDat`](@ref) constructor.
"""
function InterpolationModel(
        tbdat::WannierIO.TbDat;
        fractional_centers::AbstractVector,
        wsvec::Union{Nothing, WannierIO.WsvecDat} = nothing,
        atom_positions::AbstractVector = Vec3{Float64}[],
        atom_labels::AbstractVector = String[],
        symmetry = nothing,
    )
    return _interpolation_model_from_w90_hamiltonian(
        tbdat.H,
        tbdat.lattice,
        tbdat.Rvectors,
        tbdat.Rdegens;
        fractional_centers,
        wsvec,
        atom_positions,
        atom_labels,
        symmetry,
    )
end
