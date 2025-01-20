
function compute_phase_factor(H::Wannier.TBOperator)
    # It should be phase rotation irrespective of kpoints, nor R vectors.
    # Therefore, just look at the origin R=0.
    H0 = H[0, 0, 0]
    # L_ij := u_i^* u_j = H_ij^* / ||H_ij||
    H0_1 = H0[1, :]  # only need 1st row
    L = conj(H0_1) ./ norm.(H0_1)
    # L[1] is always 1 since H0[1,1] ∈ ℝ.
    # I assume u_1 is always 1, then set u_2, ..., u_nwann according to the
    # values of L[2:end].
    nwann = n_wannier(H)
    phases = ones(eltype(L), nwann)
    phases[2:end] = L[2:end]

    #=
    The following is optional, since one can multiply -1 to any MLWF and it
    won't change the spectrum of the Hamiltonian.
    However, to reach some sort of consistency, I enforce the sign of the
    elements of the lowest eigenvector to be real.
    The physical intuition is that we want the lowest eigenvector
    v = [v₁, ..., vₙ] = ∑ᵢ wᵢ * vᵢ (i.e., bonding orbital, having the lowest
    eigenvalue) to be constructive interference of MLWFs |wᵢ>, meaning the sign
    of the vᵢ should be positive, provided that all the |wᵢ> are consistently
    "positive" (e.g. on real space grid points, the maximum magnitude of
    equivalent MLWFs should be all positive). With this condition, the hopping
    amplitudes between equivalent MLWFs will have the same sign; otherwise,
    they can have same magnitude but different signs.
    =#
    P = Diagonal(phases)
    PHP = P' * H0 * P  # this should be real
    maximag = maximum(norm.(imag.(PHP)))
    (maximag < 1e-4) || @error "Hamiltonian cannot be real" phases maximag

    v1 = eigvecs(real(PHP))[:, 1]
    phases .*= sign.(v1)

    #=
    In the end, we have fixed the relative phase between all the MLWFs, but
    there is still a global phase rotation. However, it won't change the
    realness & sign of the Hamiltonian.
    There is no way to get this global phase, except computing it from real
    space MLWF, since it depend on the phase of Bloch functions.
    =#
    phases
end
