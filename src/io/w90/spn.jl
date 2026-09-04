export read_w90_tb_chk_spn

"""
    read_w90_tb_chk_spn(prefix; chk=prefix*".chk", spn=prefix*".spn",
                        fractional_centers=nothing, atom_positions=[],
                        atom_labels=[], symmetry=nothing)

Read Wannier90 `tb.dat`, optional `wsvec.dat`, `chk`, and `spn` data directly
into one [`InterpolationModel`](@ref). The returned model contains
`:hamiltonian`, `:berry_connection`, and `:spin` primitives on a common packed
real-space domain.
"""
function read_w90_tb_chk_spn(
        prefix::AbstractString;
        chk::AbstractString = prefix * ".chk",
        spn::AbstractString = prefix * ".spn",
        fractional_centers::Union{Nothing, AbstractVector} = nothing,
        atom_positions::AbstractVector = Vec3{Float64}[],
        atom_labels::AbstractVector = String[],
        symmetry = nothing,
    )
    tbdat = read_w90_tb_dat(prefix * "_tb.dat")
    chkdat = read_chk(chk)
    gauges = WannierIO.gauge_matrices(chkdat)
    spin_input = BlochOperator(read_spn(spn))
    _validate_bloch_operator(
        :spin, spin_input, size(gauges, 1), size(gauges, 3)
    )
    spin_kspace = _transform_bloch_operator(spin_input, gauges)
    phase_type = promote_type(eltype(spin_kspace), eltype(tbdat.H))
    phase = _quotient_fourier_phase(chkdat.kpoints, tbdat.Rvectors, phase_type)
    spin_coefficients = _quotient_fourier_coefficients(spin_kspace, phase)

    centers = isnothing(fractional_centers) ?
        _tb_fractional_centers(tbdat) : fractional_centers
    operator_coefficients = (;
        hamiltonian = tbdat.H,
        berry_connection = _tb_position(tbdat),
        spin = spin_coefficients,
    )
    descriptions = (;
        hamiltonian = (;
            law = Scalar(time_reversal = Even()),
            hermitian = true,
        ),
        berry_connection = _operator_description(BerryConnection()),
        spin = _operator_description(spin_input),
    )
    return _interpolation_model_from_w90_operators(
        operator_coefficients,
        descriptions,
        tbdat.lattice,
        tbdat.Rvectors,
        tbdat.Rdegens;
        fractional_centers = centers,
        wsvec = _read_optional_wsvec(prefix),
        atom_positions,
        atom_labels,
        symmetry,
    )
end
