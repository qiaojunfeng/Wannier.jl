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
    Hᴿ::TBHamiltonian,
    kpoints::AbstractVector{Vec3{T}};
    use_degen_pert::Bool=false,
    degen::T=1e-4,
) where {T<:Real}
    n_wann = size(Hᴿ[1].block, 1)
    n_kpts = length(kpoints)
    size(Hᴿ[1], 1) == size(Hᴿ[1], 2) || error("Hᴿ is not square")

    # first I need Hamiltonian eigenvalues and eigenvectors
    H = HamiltonianKGrid(Hᴿ, kpoints)
    E, U = H.eigvals, H.eigvecs
    # now size(E) = (n_wann, n_kpts), size(U) = (n_wann, n_wann, n_kpts)

    # velocity V = dH / dk
    # Wannier gauge V
    Vᵂ = zeros(Complex{T}, n_wann, n_wann)
    # Hamiltonian gauge V
    Vᴴ = zeros(Complex{T}, n_wann, n_wann)
    # diagonal part of Vᴴ, and is real, last index is cartesian direction
    vᴴ = zeros(T, n_wann, n_kpts, 3)

    for a in 1:3  # three Cartesian directions, a ∈ {x, y, z}
        RH = map(
            x -> TBBlock(x.R_cryst, x.R_cart, im .* x.R_cart[a] .* x.block, x.tb_block), Hᴿ
        )

        for ik in 1:n_kpts
            fill!(Vᵂ, 0)
            # for diagonal part, U† Vᵂ U = Vᴴ
            invfourier(RH, kpoints[ik]) do i, iR, Rcart, b, fac
                Vᵂ[i] += fac * b.block[i]
            end

            Uᵏ = U[ik]
            Vᴴ .= Uᵏ' * Vᵂ * Uᵏ
            vᴴ[:, ik, a] = real(diag(Vᴴ))  # YWVS Eq.27

            use_degen_pert || continue

            # now considering possible degeneracies
            mask = trues(n_wann)  # E[mask, ik] are eigenvalues to be checked

            while any(mask)
                e = E[ik][mask][1]
                # indexes of degenerate eigenvalues
                idx = abs.(E[ik] .- e) .< degen
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

"""
    get_dH_da(Rvectors, Hᴿ, kpoints)

Compute the derivative of the Hamiltonian ``H`` with respect to three Cartesian

directions.
YWVS Eq. 26.
"""
function get_dH_da(Hᴿ::TBHamiltonian, kpoints::AbstractVector{Vec3{T}}) where {T<:Real}
    E, U, Haᴴ, D = _get_D(Hᴿ, kpoints; use_degen_pert=false)
    # size(E) = (n_wann, n_kpts)
    # size(U) = (n_wann, n_wann, n_kpts)
    # size(Haᴴ) = size(D) = (n_wann, n_wann, n_kpts, 3)

    # also including non-diagonal part
    n_kpts = length(kpoints)
    n_wann = length(E[1])
    dH_da = zeros(Complex{T}, n_wann, n_wann, n_kpts, 3)

    HDa = zeros(Complex{T}, n_wann, n_wann)

    for a in 1:3
        for ik in 1:n_kpts
            HDa .= Diagonal(E[ik]) * D[:, :, ik, a]
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
function get_d2H_dadb(Hᴿ::TBHamiltonian, kpoints::AbstractVector{Vec3{T}}) where {T<:Real}
    n_wann = size(Hᴿ[1], 1)
    n_kpts = length(kpoints)
    size(Hᴿ[1], 1) == size(Hᴿ[1], 2) || error("Hᴿ is not square")

    E, U, Haᴴ, D = _get_D(Hᴿ, kpoints; use_degen_pert=false)
    # size(E) = (n_wann, n_kpts)
    # size(U) = (n_wann, n_wann, n_kpts)
    # size(Haᴴ) = size(D) = (n_wann, n_wann, n_kpts, 3)

    # 2nd order derivative Hab = dv / dk = d²H / dk²
    # Wannier gauge Hab, last two indexes are Cartesian directions
    Habᵂ = zeros(Complex{T}, n_wann, n_wann)
    # Hamiltonian gauge Hab, actually the covariant part of Hab
    # i.e., the ``\bar{H}_{\alpha\beta}^{(H)}`` in YWVS Eq. 28
    Habᴴ = zeros(Complex{T}, n_wann, n_wann, n_kpts, 3, 3)

    # to cartesian in angstrom
    for a in 1:3  # three Cartesian directions, a ∈ {x, y, z}
        for b in 1:3
            RaRbH = map(
                x -> TBBlock(
                    x.R_cryst,
                    x.R_cart,
                    -x.R_cart[a] * x.R_cart[b] .* x.block,
                    x.tb_block,
                ),
                Hᴿ,
            )

            for ik in 1:n_kpts
                fill!(Habᵂ, 0)

                invfourier(RaRbH, kpoints[ik]) do i, iR, Rcart, b, fac
                    Habᵂ[i] += fac * b.block[i]
                end

                Uᵏ = U[ik]
                Habᴴ[:, :, ik, a, b] .= Uᵏ' * Habᵂ * Uᵏ

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
    Rvectors::RVectorsMDRS{T}, Hᴿ::TBHamiltonian, kpoints::AbstractVector{Vec3{T}}; dk=1e-3
) where {T<:Real}
    n_wann = size(Hᴿ[1], 1)
    n_kpts = length(kpoints)
    size(Hᴿ[1], 1) == size(Hᴿ[1], 2) || error("Hᴿ is not square")

    # the final velocity along 3 Cartesian directions
    V = zeros(T, n_wann, n_kpts, 3)
    # the 6 interpolated kpoints, each column is a fractional coordinates
    Δk = [
        Vec3(-dk, 0, 0),
        Vec3(dk, 0, 0),
        Vec3(0, -dk, 0),
        Vec3(0, dk, 0),
        Vec3(0, 0, -dk),
        Vec3(0, 0, dk),
    ]

    # the interpolated Hamiltonain
    Hᵏ = zeros(Complex{T}, n_wann, n_wann)
    # the interpolated eigenvalues
    E = zeros(T, n_wann, 6)
    # to fractional
    recip_latt = reciprocal_lattice(Rvectors.lattice)
    Δkᶠ = map(k -> inv(recip_latt) * k, Δk)

    for ik in 1:n_kpts
        for id in 1:6
            fill!(Hᵏ, 0)

            invfourier(Hᴿ, kpoints[ik] .+ Δkᶠ[id]) do i, iR, Rcart, b, fac
                Hᵏ[i] += fac * b.block[i]
            end

            E[:, id] .= real(eigen(Hᵏ).values)  # only eigenvalues are needed
        end

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
    Rvectors::RVectorsMDRS{T}, Hᴿ::TBHamiltonian, kpoints::AbstractVector{Vec3{T}}; dk=1e-3
) where {T<:Real}
    n_wann = size(Hᴿ[1].block, 1)
    n_kpts = length(kpoints)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    size(Hᴿ[1].block, 1) == size(Hᴿ[1].block, 2) || error("Hᴿ is not square")
    length(Hᴿ) == n_r̃vecs || error("Hᴿ has wrong size")

    # the final velocity along 3 x 3 Cartesian directions
    μ = zeros(T, n_wann, n_kpts, 3, 3)
    # the 1 + 6 + 12 interpolated kpoints, each column is a fractional coordinates
    # note the original kpoint is included
    # 6 points at k ± dk, for computing dϵ/dk at k ± dk/2, along the same direction as dk
    # 12 points at k ± dk/2, for computing dϵ/dk at k ± dk/2, but along other 2 directions
    Δk_ = [
        0 -dk dk 0 0 0 0 -dk/2 -dk/2 -dk/2 -dk/2 dk/2 dk/2 dk/2 dk/2 0 0 0 0
        0 0 0 -dk dk 0 0 -dk/2 dk/2 0 0 -dk/2 dk/2 0 0 -dk/2 -dk/2 dk/2 dk/2
        0 0 0 0 0 -dk dk 0 0 -dk/2 dk/2 0 0 -dk/2 dk/2 -dk/2 dk/2 -dk/2 dk/2
    ]
    Δk = [Vec3(Δk_[:, ik]) for ik in 1:size(Δk_, 2)]
    # the interpolated Hamiltonain
    Hᵏ = zeros(Complex{T}, n_wann, n_wann)
    # the interpolated eigenvalues
    E = [zeros(T, n_wann) for i in 1:19]
    # dϵ/dk at 6 points: k ± dk/2, along 3 Cartesian directions
    dEdk = zeros(T, n_wann, 6, 3)
    # to fractional
    recip_latt = reciprocal_lattice(Rvectors.lattice)
    Δkᶠ = map(k -> inv(recip_latt) * k, Δk)
    for ik in 1:n_kpts
        for i2 in 1:19
            fill!(Hᵏ, 0)
            invfourier(Hᴿ, kpoints[ik] .+ Δkᶠ[i2]) do i, iR, Rcart, b, fac
                Hᵏ[i] += fac * b.block[i]
            end
            E[i2] .= real(eigen(Hᵏ).values)  # only eigenvalues are needed
        end
        # dϵ/dk at k ± dk/2, along the same direction as dk
        dEdk[:, 1, 1] .= (E[2] - E[1]) / -dk
        dEdk[:, 2, 1] .= (E[3] - E[1]) / dk
        dEdk[:, 3, 2] .= (E[4] - E[1]) / -dk
        dEdk[:, 4, 2] .= (E[5] - E[1]) / dk
        dEdk[:, 5, 3] .= (E[6] - E[1]) / -dk
        dEdk[:, 6, 3] .= (E[7] - E[1]) / dk

        # dϵ/dk at k ± dk/2, but along other 2 directions
        dEdk[:, 1, 2] .= (E[9] - E[8]) / dk
        dEdk[:, 1, 3] .= (E[11] - E[10]) / dk
        dEdk[:, 2, 2] .= (E[13] - E[12]) / dk
        dEdk[:, 2, 3] .= (E[15] - E[14]) / dk
        dEdk[:, 3, 1] .= (E[12] - E[8]) / dk
        dEdk[:, 3, 3] .= (E[17] - E[16]) / dk
        dEdk[:, 4, 1] .= (E[13] - E[9]) / dk
        dEdk[:, 4, 3] .= (E[19] - E[18]) / dk
        dEdk[:, 5, 1] .= (E[14] - E[10]) / dk
        dEdk[:, 5, 2] .= (E[18] - E[16]) / dk
        dEdk[:, 6, 1] .= (E[15] - E[11]) / dk
        dEdk[:, 6, 2] .= (E[19] - E[17]) / dk

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
