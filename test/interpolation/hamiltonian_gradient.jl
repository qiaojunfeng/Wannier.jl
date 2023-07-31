@testitem "VelocityInterpolator FiniteDifferenceVelocity" begin
    using Wannier.Datasets
    hamiltonian, _ = read_w90_tb(dataset"Si2_valence/reference/MDRS/Si2_valence")
    interp = Wannier.VelocityInterpolator(hamiltonian)

    # choose a random k-point such that there is no degeneracy
    k = [0.1, 0.2, 0.3]
    V = interp(k, Wannier.FiniteDifferenceVelocity())

    ref_V = [
        [2.9658709201356714, 1.3924552106434618, -1.8952333036281743e-6],
        [-5.823622720557475, -4.4448034497241995, 1.4344967436130673e-6],
        [-4.132392387731443, -4.20084007884336, 1.2558989403999021e-6],
        [-4.956384133104397, 2.6956732295380093, 1.9835657560918207e-6],
    ]
    @test all(isapprox.(V, ref_V; atol=1e-7))
end

@testitem "VelocityInterpolator FourierSpaceVelocity" begin
    using Wannier.Datasets
    hamiltonian, _ = read_w90_tb(dataset"Si2_valence/reference/MDRS/Si2_valence")
    interp = Wannier.VelocityInterpolator(hamiltonian)

    k = [0.1, 0.2, 0.3]
    V = interp(k, Wannier.FourierSpaceVelocity())

    ref_V = interp(k, Wannier.FiniteDifferenceVelocity())
    @test all(isapprox.(V, ref_V; atol=5e-5))
end

@testitem "HamiltonianGradientInterpolator" begin
    using LinearAlgebra
    using Wannier.Datasets
    hamiltonian, _ = read_w90_tb(dataset"Si2_valence/reference/MDRS/Si2_valence")
    interp = Wannier.HamiltonianGradientInterpolator(hamiltonian)

    k = [0.1, 0.2, 0.3]
    dH = interp(k)

    ref_dHdiag = Wannier.VelocityInterpolator(hamiltonian)(
        k, Wannier.FourierSpaceVelocity()
    )
    @test diag(dH) ≈ ref_dHdiag
    V_offdiag = dH - diagm(ref_dHdiag)
    # element of dH are complex MVec3, so this also ensure that
    # the imaginary part is zero
    @test all(isapprox.(V_offdiag, Ref([0, 0, 0]); atol=1e-10))
end
