using LinearAlgebra

"""
Compute Wilson loop along k_y for a 2D Hamiltonian (assume vaccume along z).
"""
# function wilson_loop(hamiltonian::TBOperator, centers::AbstractVector)
    nkx = 20
    nky = 20
    nocc = 10

    cd("/Users/junfeng/git/wannier_tools/examples/Bi2Se3")
    wout = read_wout("wannier90.wout")
    hamiltonian = read_w90_hr("wannier90", wout.lattice)
    interp = HamiltonianInterpolator(hamiltonian)
    centers = wout.centers

    nwann = length(centers)
    kpts_x = range(start=0, step=1/nkx, length=nkx)
    Δky = 1/nky

    loops = map(enumerate(kpts_x)) do ikx, kx
        kpoints = [Vec3(kx, ky, 0.0) for ky in range(start=0, step=Δky, length=nky)]
        E, V = interp(kpoints)
        # check if occupied bands are isolated from unoccupied bands
        Egap = [e[occ+1] - e[nocc] for e in E]
        idx_gap = argmin(Egap)
        gap = Egap[idx_gap]
        if gap < 0
            @warn "Energy gap $gap < 0 at (kx, ky, kz) = $(kpoints[idx_gap])"
        end
        # The Wannier-gauge overlap matrix Mᵂₖ₁ₖ₂ = <uᵂₖ₁|uᵂₖ₂>
        # = <ψᵂₖ₁|exp(-iΔk⋅r̂)|ψᵂₖ₂>  where Δk = k₂-k₁, r̂ is the position operator
        # = ∑_R₁ exp(-ik₁R₁) <WF(R₁)| * exp(-iΔk⋅r̂) * ∑_R₂ exp(ik₂R₂) |WF(R₂)>
        # = ∑_R₁∑_R₂ exp(i(k₂R₂-k₁R₁)) <WF(R₁)|exp(-iΔk⋅r̂)|WF(R₂)>
        # = ∑_R₁∑_R₂ exp(i(k₂R₂-k₁R₁)) <WF(R₁)|exp(-iΔk⋅(r̄+R₂))|WF(R₂)>  where r̄ is the center of |WF(R₀)>
        # = ∑_R₁∑_R₂ exp(i(k₂R₂-k₁R₁)) exp(-iΔk⋅(r̄+R₂)) <WF(R₁)|WF(R₂)>
        # = ∑_R exp(i Δk R) exp(-iΔk⋅(r̄+R))  due to orthogonality of WF
        # = N exp(-iΔk⋅r̄)
        Mᵂ = diag(exp(-im * [0.0, Δky, 0.0] .⋅ centers))
        # to Bloch gauge
        Mᴮ = Vocc' * Mᵂ * Vocc
        # the Wilson loop at this kx line
        loop_kx = map(1:nky) do iky
            Vocc = V[iky][: 1:nocc]
    end
    loop
# end
