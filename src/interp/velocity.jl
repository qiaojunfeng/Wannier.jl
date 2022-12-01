using LinearAlgebra: diag

"""
    velocity(Rvectors, Hᴿ, kpoints; use_degen_pert=false, degen=1e-4)

Compute velocity along three Cartesian directions.

# Arguments
- `Rvectors`: `RVectorsMDRS`
- `Hᴿ`: `n_wann * n_wann * n_r̃vecs`, the Hamiltonian matrix
- `kpoints`: `3 * n_kpts`, each column is a fractional kpoints coordinates

# Keyword arguments
- `use_degen_pert`: use perturbation treatment for degenerate eigenvalues
- `degen`: degeneracy threshold in eV

# Return
- `E`: `n_wann * n_kpts`, energy eigenvalues
- `v`: `n_wann * n_kpts * 3`, velocity along three Cartesian directions

!!! warning

    `Wannier90` by default set `use_degen_pert = false`.
    Switching off perturbation treatment for degenerate eigenvalues will
    lead to somewhat arbitrary velocities for degenerate states.
    Thus I set `use_degen_pert = true` by default.
"""
function velocity(
    Rvectors::RVectorsMDRS{T},
    Hᴿ::Array{Complex{T},3},
    kpoints::AbstractMatrix{T};
    use_degen_pert::Bool=true,
    degen::T=1e-4,
) where {T<:Real}
    n_wann = size(Hᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = Rvectors.n_r̃vecs  # only MDRSv2 is implemented for the moment
    @assert size(Hᴿ, 1) == size(Hᴿ, 2)
    @assert size(Hᴿ, 3) == n_r̃vecs

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

    for iv in 1:3
        Riv = reshape(Rᶜ[iv, :], 1, 1, n_r̃vecs)
        RH = im * Riv .* Hᴿ
        Vᵂ .= invfourier(Rvectors, RH, kpoints)

        for ik in 1:n_kpts
            # for diagonal part, U† Vᵂ U = Vᴴ
            Uᵏ = @view U[:, :, ik]
            Vᴴ .= Uᵏ' * Vᵂ[:, :, ik] * Uᵏ
            vᴴ[:, ik, iv] = real(diag(Vᴴ))  # YWVS Eq.27

            use_degen_pert || continue

            # now considering possible degeneracies
            mask = trues(n_wann)  # E[mask, ik] are eigenvalues to be checked
            while any(mask)
                e = E[mask, ik][1]
                # indexes of degenerate eigenvalues
                idx = abs.(E[:, ik] .- e) .< degen
                if count(idx) > 1
                    # diagonalize the submatrix
                    vᴴ[idx, ik, iv] .= real(eigen(Vᴴ[idx, idx]).values)
                end

                mask[idx] .= false
            end
        end
    end

    return E, vᴴ
end
