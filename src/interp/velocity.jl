using LinearAlgebra

"""
    velocity(Rvectors, Hᴿ, kpoints; use_degen_pert=false, degen=1e-4)

Compute velocity along three Cartesian directions.

YWVS Eq. 27.

# Arguments
- `Rvectors`: `RVectorsMDRS`
- `Hᴿ`: `n_wann * n_wann * n_r̃vecs`, the Hamiltonian matrix
- `kpoints`: `3 * n_kpts`, each column is a fractional kpoints coordinates

# Keyword arguments
- `use_degen_pert`: use perturbation treatment for degenerate eigenvalues
- `degen`: degeneracy threshold in eV

# Return
- `E`: `n_wann * n_kpts`, energy eigenvalues
- `v`: `n_wann * n_kpts * 3`, velocity along three Cartesian directions,
    in unit `hbar * m / s`

!!! warning

    `Wannier90` by default set `use_degen_pert = false`.
    In 3D, and for N degenerate states, the velocity is a tensor
    of size N * N * 3, where 3 is for three Cartesian directions.
    Thus I cannot simultaneously diagonalize the tensor for all 3 directions.
    This means I can only use perturbation treatment for one of the directions,
    and only in that direction the velocity matrix is diagonal.
    So for degenerate states, the velocity is not well defined, and the results
    are meaningless, instead one should use the full velocity matrix which
    also include non-diagonal part, see [`get_dH_da`](@ref).
"""
function velocity(
    Rvectors::RVectorsMDRS{T},
    Hᴿ::Array{Complex{T},3},
    kpoints::AbstractMatrix{T};
    use_degen_pert::Bool=false,
    degen::T=1e-4,
) where {T<:Real}
    n_wann = size(Hᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    size(Hᴿ, 1) == size(Hᴿ, 2) || error("Hᴿ is not square")
    size(Hᴿ, 3) == n_r̃vecs || error("Hᴿ has wrong size")

    # first I need Hamiltonian eigenvalues and eigenvectors
    H = invfourier(Rvectors, Hᴿ, kpoints)
    E, U = diag_Hk(H)
    # now size(E) = (n_wann, n_kpts), size(U) = (n_wann, n_wann, n_kpts)

    # velocity V = dH / dk
    # Wannier gauge V
    Vᵂ = zeros(Complex{T}, n_wann, n_wann, n_kpts)
    # Hamiltonian gauge V
    Vᴴ = zeros(Complex{T}, n_wann, n_wann)
    # diagonal part of Vᴴ, and is real, last index is cartesian direction
    vᴴ = zeros(T, n_wann, n_kpts, 3)

    # to cartesian in angstrom
    Rᶜ = Rvectors.lattice * Rvectors.R̃vectors.R

    for a in 1:3  # three Cartesian directions, a ∈ {x, y, z}
        Ra = reshape(Rᶜ[a, :], 1, 1, n_r̃vecs)
        RH = im * Ra .* Hᴿ
        invfourier!(Vᵂ, Rvectors, RH, kpoints)

        for ik in 1:n_kpts
            # for diagonal part, U† Vᵂ U = Vᴴ
            Uᵏ = @view U[:, :, ik]
            Vᴴ .= Uᵏ' * Vᵂ[:, :, ik] * Uᵏ
            vᴴ[:, ik, a] = real(diag(Vᴴ))  # YWVS Eq.27

            use_degen_pert || continue

            # now considering possible degeneracies
            mask = trues(n_wann)  # E[mask, ik] are eigenvalues to be checked
            while any(mask)
                e = E[mask, ik][1]
                # indexes of degenerate eigenvalues
                idx = abs.(E[:, ik] .- e) .< degen
                if count(idx) > 1
                    # diagonalize the submatrix
                    vᴴ[idx, ik, a] .= real(eigen(Vᴴ[idx, idx]).values)
                end

                mask[idx] .= false
            end
        end
    end

    return E, vᴴ
end

@doc raw"""
    _get_D(Rvectors, Hᴿ, kpoints; use_degen_pert=false, degen=1e-4)

Compute the matrix D in YWVS Eq. 25 (or Eq. 32 if `use_degen_pert = true`).

# Arguments
- `Rvectors`: `RVectorsMDRS`
- `Hᴿ`: `n_wann * n_wann * n_r̃vecs`, the Hamiltonian matrix
- `kpoints`: `3 * n_kpts`, each column is a fractional kpoints coordinates

# Keyword arguments
- `use_degen_pert`: use perturbation treatment for degenerate eigenvalues
- `degen`: degeneracy threshold in eV

# Return
- `E`: `n_wann * n_kpts`, energy eigenvalues
- `U`: `n_wann * n_wann * n_kpts`, the unitary transformation matrix
- `Haᴴ`: `n_wann * n_wann * n_kpts * 3`, the covariant part of derivative of
    Hamiltonian in Hamiltonian gauge,
    the ``\bar{H}_{\alpha}^{(H)}`` in YWVS Eq. 26
- `D`: `n_wann * n_wann * n_kpts * 3`, the matrix D in YWVS Eq. 25 or Eq. 32

!!! warning

    If `use_degen_pert = true`, the degenerate subspace is rotated such that
    ``\bar{H}_{x}^{(H)}`` is diagonal, note only the ``x`` direction.
    I cannot diagonalize simultaneously all the three directions.
"""
function _get_D(
    Rvectors::RVectorsMDRS{T},
    Hᴿ::Array{Complex{T},3},
    kpoints::AbstractMatrix{T};
    use_degen_pert::Bool=false,
    degen::T=1e-4,
) where {T<:Real}
    n_wann = size(Hᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    size(Hᴿ, 1) == size(Hᴿ, 2) || error("Hᴿ is not square")
    size(Hᴿ, 3) == n_r̃vecs || error("Hᴿ has wrong size")

    # first I need Hamiltonian eigenvalues and eigenvectors
    H = invfourier(Rvectors, Hᴿ, kpoints)
    E, U = diag_Hk(H)
    # now size(E) = (n_wann, n_kpts), size(U) = (n_wann, n_wann, n_kpts)

    # velocity Ha = dH / dk_a
    # Wannier gauge Ha
    Haᵂ = zeros(Complex{T}, n_wann, n_wann, n_kpts)
    # the covariant part of Hamiltonian gauge Ha, i.e., Haᴴ = U† Haᵂ U
    # also the ``\bar{H}_{\alpha}^{(H)}`` in YWVS Eq. 26
    # last index for the three Cartesian directions
    Haᴴ = zeros(Complex{T}, n_wann, n_wann, n_kpts, 3)
    # The D matrix = U† ∂U in Hamiltonian gauge, i.e. YWVS Eq. 25 or Eq. 32
    D = zeros(Complex{T}, n_wann, n_wann, n_kpts, 3)

    # to cartesian in angstrom
    Rᶜ = Rvectors.lattice * Rvectors.R̃vectors.R

    for a in 1:3  # three Cartesian directions, a ∈ {x, y, z}
        Ra = reshape(Rᶜ[a, :], 1, 1, n_r̃vecs)
        RH = im * Ra .* Hᴿ
        invfourier!(Haᵂ, Rvectors, RH, kpoints)

        for ik in 1:n_kpts
            # Haᴴ = U† Haᵂ U
            Uᵏ = @view U[:, :, ik]
            Haᴴ[:, :, ik, a] .= Uᵏ' * Haᵂ[:, :, ik] * Uᵏ

            # the D matrix
            ΔE = E[:, ik] .- E[:, ik]'
            # assign a nonzero number to the diagonal elements for inversion
            ΔE[diagind(ΔE)] .= 1
            Dka = @view D[:, :, ik, a]
            Dka .= Haᴴ[:, :, ik, a] ./ (-ΔE)
            Dka[diagind(Dka)] .= 0

            # TODO: maybe it is helpful to run at least once the perturbation treatment
            # for one Cartesian direction, to avoid vanishing denominator in D matrix
            use_degen_pert || continue

            # now considering possible degeneracies
            mask = trues(n_wann)  # E[mask, ik] are eigenvalues to be checked
            while any(mask)
                e = E[mask, ik][1]
                # indexes of degenerate eigenvalues
                idx = abs.(E[:, ik] .- e) .< degen
                if count(idx) > 1
                    # I can only run once the diagonalization for only one Cartesian
                    # direction, and update the U matrix. The following directions
                    # will use the updated U matrix, and I only set the D matrix to
                    # zero for the degenerate subspace.
                    if a == 1
                        # diagonalize the submatrix
                        v, u = eigen(Haᴴ[idx, idx, ik, a])
                        Haᴴ[idx, idx, ik, a] .= 0
                        diag(Haᴴ[idx, idx, ik, a]) .= v
                        # update U such that in Hamiltonain gauge both H and Ha
                        # are diagonal in the degenerate subspace
                        U[idx, idx, ik] *= u
                    end
                    # the D matrix
                    D[idx, idx, ik, a] .= 0
                end

                mask[idx] .= false
            end
        end
    end

    return E, U, Haᴴ, D
end

"""
    get_dH_da(Rvectors, Hᴿ, kpoints)

Compute the derivative of the Hamiltonian ``H`` with respect to three Cartesian
directions.

YWVS Eq. 26.
"""
function get_dH_da(
    Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}, kpoints::AbstractMatrix{T}
) where {T<:Real}
    n_wann = size(Hᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    size(Hᴿ, 1) == size(Hᴿ, 2) || error("Hᴿ is not square")
    size(Hᴿ, 3) == n_r̃vecs || error("Hᴿ has wrong size")

    E, U, Haᴴ, D = _get_D(Rvectors, Hᴿ, kpoints; use_degen_pert=false)
    # size(E) = (n_wann, n_kpts)
    # size(U) = (n_wann, n_wann, n_kpts)
    # size(Haᴴ) = size(D) = (n_wann, n_wann, n_kpts, 3)

    # also including non-diagonal part
    dH_da = zeros(Complex{T}, n_wann, n_wann, n_kpts, 3)

    HDa = zeros(Complex{T}, n_wann, n_wann)

    for a in 1:3
        for ik in 1:n_kpts
            HDa .= Diagonal(E[:, ik]) * D[:, :, ik, a]
            dH_da[:, :, ik, a] = Haᴴ[:, :, ik, a] + HDa + HDa'
        end
    end

    return dH_da
end

"""
    get_d2H_dadb(Rvectors, Hᴿ, kpoints)

Compute the second derivative of the Hamiltonian ``H`` with respect to three
Cartesian directions.

YWVS Eq. 28.
"""
function get_d2H_dadb(
    Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}, kpoints::AbstractMatrix{T}
) where {T<:Real}
    n_wann = size(Hᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    size(Hᴿ, 1) == size(Hᴿ, 2) || error("Hᴿ is not square")
    size(Hᴿ, 3) == n_r̃vecs || error("Hᴿ has wrong size")

    E, U, Haᴴ, D = _get_D(Rvectors, Hᴿ, kpoints; use_degen_pert=false)
    # size(E) = (n_wann, n_kpts)
    # size(U) = (n_wann, n_wann, n_kpts)
    # size(Haᴴ) = size(D) = (n_wann, n_wann, n_kpts, 3)

    # 2nd order derivative Hab = dv / dk = d²H / dk²
    # Wannier gauge Hab, last two indexes are Cartesian directions
    Habᵂ = zeros(Complex{T}, n_wann, n_wann, n_kpts)
    # Hamiltonian gauge Hab, actually the covariant part of Hab
    # i.e., the ``\bar{H}_{\alpha\beta}^{(H)}`` in YWVS Eq. 28
    Habᴴ = zeros(Complex{T}, n_wann, n_wann, n_kpts, 3, 3)

    # to cartesian in angstrom
    Rᶜ = Rvectors.lattice * Rvectors.R̃vectors.R

    for a in 1:3  # three Cartesian directions, a ∈ {x, y, z}
        Ra = reshape(Rᶜ[a, :], 1, 1, n_r̃vecs)
        for b in 1:3
            Rb = reshape(Rᶜ[b, :], 1, 1, n_r̃vecs)
            RaRbH = -Ra .* Rb .* Hᴿ  # YWVS Eq. 30
            invfourier!(Habᵂ, Rvectors, RaRbH, kpoints)

            for ik in 1:n_kpts
                Uᵏ = @view U[:, :, ik]
                Habᴴ[:, :, ik, a, b] .= Uᵏ' * Habᵂ[:, :, ik] * Uᵏ

                HaDb = Haᴴ[:, :, ik, a] * D[:, :, ik, b]
                Habᴴ[:, :, ik, a, b] += HaDb + HaDb'  # YWVS Eq. 28
            end
        end
    end

    return Habᴴ
end

"""
    velocity_fd(Rvectors, Hᴿ, kpoints; Δk=1e-3)

Compute the velocity using finite differences of 2nd order.

PRB 93, 205147 (2016)  Eq. 80.
"""
function velocity_fd(
    Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}, kpoints::AbstractMatrix{T}; dk=1e-3
) where {T<:Real}
    n_wann = size(Hᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    size(Hᴿ, 1) == size(Hᴿ, 2) || error("Hᴿ is not square")
    size(Hᴿ, 3) == n_r̃vecs || error("Hᴿ has wrong size")

    # the final velocity along 3 Cartesian directions
    V = zeros(T, n_wann, n_kpts, 3)
    # the 6 interpolated kpoints, each column is a fractional coordinates
    Δk = [
        -dk dk 0 0 0 0
        0 0 -dk dk 0 0
        0 0 0 0 -dk dk
    ]
    # the interpolated Hamiltonain
    Hᵏ = zeros(Complex{T}, n_wann, n_wann, 6)
    # the interpolated eigenvalues
    E = zeros(T, n_wann, 6)
    # to fractional
    recip_latt = get_recip_lattice(Rvectors.lattice)
    Δkᶠ = inv(recip_latt) * Δk

    for ik in 1:n_kpts
        invfourier!(Hᵏ, Rvectors, Hᴿ, kpoints[:, ik] .+ Δkᶠ)
        E .= diag_Hk(Hᵏ)[1]  # only eigenvalues are needed

        V[:, ik, 1] = (E[:, 2] - E[:, 1]) / (2dk)
        V[:, ik, 2] = (E[:, 4] - E[:, 3]) / (2dk)
        V[:, ik, 3] = (E[:, 6] - E[:, 5]) / (2dk)
    end

    return V
end

"""
    effmass_fd(Rvectors, Hᴿ, kpoints; Δk=1e-3)

Compute the inverse of effective mass using finite differences of 2nd order.

Apply twice PRB 93, 205147 (2016)  Eq. 80.
"""
function effmass_fd(
    Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}, kpoints::AbstractMatrix{T}; dk=1e-3
) where {T<:Real}
    n_wann = size(Hᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    size(Hᴿ, 1) == size(Hᴿ, 2) || error("Hᴿ is not square")
    size(Hᴿ, 3) == n_r̃vecs || error("Hᴿ has wrong size")

    # the final velocity along 3 x 3 Cartesian directions
    μ = zeros(T, n_wann, n_kpts, 3, 3)
    # the 1 + 6 + 12 interpolated kpoints, each column is a fractional coordinates
    # note the original kpoint is included
    # 6 points at k ± dk, for computing dϵ/dk at k ± dk/2, along the same direction as dk
    # 12 points at k ± dk/2, for computing dϵ/dk at k ± dk/2, but along other 2 directions
    Δk = [
        0 -dk dk 0 0 0 0 -dk/2 -dk/2 -dk/2 -dk/2 dk/2 dk/2 dk/2 dk/2 0 0 0 0
        0 0 0 -dk dk 0 0 -dk/2 dk/2 0 0 -dk/2 dk/2 0 0 -dk/2 -dk/2 dk/2 dk/2
        0 0 0 0 0 -dk dk 0 0 -dk/2 dk/2 0 0 -dk/2 dk/2 -dk/2 dk/2 -dk/2 dk/2
    ]
    # the interpolated Hamiltonain
    Hᵏ = zeros(Complex{T}, n_wann, n_wann, 19)
    # the interpolated eigenvalues
    E = zeros(T, n_wann, 19)
    # dϵ/dk at 6 points: k ± dk/2, along 3 Cartesian directions
    dEdk = zeros(T, n_wann, 6, 3)
    # to fractional
    recip_latt = get_recip_lattice(Rvectors.lattice)
    Δkᶠ = inv(recip_latt) * Δk

    for ik in 1:n_kpts
        invfourier!(Hᵏ, Rvectors, Hᴿ, kpoints[:, ik] .+ Δkᶠ)
        E .= diag_Hk(Hᵏ)[1]  # only eigenvalues are needed

        # dϵ/dk at k ± dk/2, along the same direction as dk
        dEdk[:, 1, 1] .= (E[:, 2] - E[:, 1]) / -dk
        dEdk[:, 2, 1] .= (E[:, 3] - E[:, 1]) / dk
        dEdk[:, 3, 2] .= (E[:, 4] - E[:, 1]) / -dk
        dEdk[:, 4, 2] .= (E[:, 5] - E[:, 1]) / dk
        dEdk[:, 5, 3] .= (E[:, 6] - E[:, 1]) / -dk
        dEdk[:, 6, 3] .= (E[:, 7] - E[:, 1]) / dk

        # dϵ/dk at k ± dk/2, but along other 2 directions
        dEdk[:, 1, 2] .= (E[:, 9] - E[:, 8]) / dk
        dEdk[:, 1, 3] .= (E[:, 11] - E[:, 10]) / dk
        dEdk[:, 2, 2] .= (E[:, 13] - E[:, 12]) / dk
        dEdk[:, 2, 3] .= (E[:, 15] - E[:, 14]) / dk
        dEdk[:, 3, 1] .= (E[:, 12] - E[:, 8]) / dk
        dEdk[:, 3, 3] .= (E[:, 17] - E[:, 16]) / dk
        dEdk[:, 4, 1] .= (E[:, 13] - E[:, 9]) / dk
        dEdk[:, 4, 3] .= (E[:, 19] - E[:, 18]) / dk
        dEdk[:, 5, 1] .= (E[:, 14] - E[:, 10]) / dk
        dEdk[:, 5, 2] .= (E[:, 18] - E[:, 16]) / dk
        dEdk[:, 6, 1] .= (E[:, 15] - E[:, 11]) / dk
        dEdk[:, 6, 2] .= (E[:, 19] - E[:, 17]) / dk

        # d²ϵ/dk² at k
        μ[:, ik, 1, 1] .= (dEdk[:, 2, 1] - dEdk[:, 1, 1]) / dk
        μ[:, ik, 1, 2] .= (dEdk[:, 4, 1] - dEdk[:, 3, 1]) / dk
        μ[:, ik, 1, 3] .= (dEdk[:, 6, 1] - dEdk[:, 5, 1]) / dk
        μ[:, ik, 2, 1] .= (dEdk[:, 2, 2] - dEdk[:, 1, 2]) / dk
        μ[:, ik, 2, 2] .= (dEdk[:, 4, 2] - dEdk[:, 3, 2]) / dk
        μ[:, ik, 2, 3] .= (dEdk[:, 6, 2] - dEdk[:, 5, 2]) / dk
        μ[:, ik, 3, 1] .= (dEdk[:, 2, 3] - dEdk[:, 1, 3]) / dk
        μ[:, ik, 3, 2] .= (dEdk[:, 4, 3] - dEdk[:, 3, 3]) / dk
        μ[:, ik, 3, 3] .= (dEdk[:, 6, 3] - dEdk[:, 5, 3]) / dk
    end

    return μ
end
