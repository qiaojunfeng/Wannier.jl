
"""
    fermi_dirac(E, ϵF, T)

Fermi-Dirac distribution function.

# Arguments
- `E`: energy eigenvalues, in eV
- `ϵF`: Fermi energy, in eV
- `T`: temperature, in Kelvin
"""
function fermi_dirac(E::AbstractVector{I}, ϵF::Real, T::Real) where {I<:Real}
    # TODO 0/0 returns NaN
    return 1 ./ (1 .+ exp.((E .- ϵF) ./ (kB * T)))
end

"""
    occupation(E, ϵF, T)

Occupations for eigenvalues.

# Arguments
- `E`: energy eigenvalues, in eV, size `(n_bands, n_kpts)`
- `ϵF`: Fermi energy, in eV
- `T`: temperature, in Kelvin

# Returns
- `occ`: occupations, size `(n_bands, n_kpts)`
"""
function occupation(E::AbstractMatrix{I}; ϵF::Real, T::Real) where {I<:Real}
    n_bands, n_kpts = size(E)
    occ = zeros(eltype(E), n_bands, n_kpts)
    for (ik, Ek) in enumerate(eachcol(E))
        occ[:, ik] = fermi_dirac(Ek, ϵF, T)
    end
    return occ
end

@doc raw"""
    occupation_matrix(model, ϵF, T)

Occupation matrix of WFs.

```math
\langle w_{n \mathbf{0}} | \hat{P} | w_{m \mathbf{0}} \rangle
```
where ``\hat{P}`` is the occupation operator.
"""
function occupation_matrix(model::Model, ϵF::Real, T::Real)
    occ = occupation(model.E; ϵF, T)
    P = get_Hk(occ, model.U)
    return dropdims(sum(P; dims=3); dims=3) / model.n_kpts
end

@doc raw"""
    occupation_matrix(model, ϵF, T)

Occupation matrix in R-space.

```math
\langle w_{n \mathbf{0}} | \hat{P} | w_{m \mathbf{R}} \rangle
```
where ``\hat{P}`` is the occupation operator.
"""
function occupation_matrix(model::InterpModel, P::AbstractArray3)
    Pᴿ = fourier(model.kRvectors, P)
    return Pᴿ
end

@doc raw"""
    charge(model, ϵF, T)

Charge of WFs.

For a WF ``w_{n \mathbf{0}}``, its charge is
```math
\langle w_{n \mathbf{0}} | \hat{P} | w_{n \mathbf{0}} \rangle
```
where ``\hat{P}`` is the occupation operator.
"""
function charge(model::Model, ϵF::Real, T::Real)
    P = occupation_matrix(model, ϵF, T)
    return real(diag(P))
end

@doc raw"""
    magnetic_moment(model, ϵF, T)

Magnetic moment of each WF.

# Arguments
- `model`: a `MagModel` in COWF or COWF+C gauge
- `ϵF`: Fermi energy, in eV
- `T`: temperature, in Kelvin

!!! warning

    This function returns sensible results only for COWF or COWF+C gauge.
"""
function magnetic_moment(model::MagModel, ϵF::Real, T::Real)
    Q_up = charge(model.up, ϵF, T)
    Q_dn = charge(model.dn, ϵF, T)
    return Q_up - Q_dn
end
