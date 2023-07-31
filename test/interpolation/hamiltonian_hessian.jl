@testitem "EffectiveMassInterpolator FiniteDifferenceEffectiveMass" begin
    using Wannier.Datasets
    hamiltonian, _ = read_w90_tb(dataset"Si2_valence/reference/MDRS/Si2_valence")
    interp = Wannier.EffectiveMassInterpolator(hamiltonian)

    # choose a random k-point such that there is no degeneracy
    k = [0.1, 0.2, 0.3]
    μ = interp(k, Wannier.FiniteDifferenceEffectiveMass())

    ref_μ = [
        [
            5.866671990872874 -0.9471753461554044 5.830003146911622e-6
            -0.9471753461554044 6.128688220030654 -1.5685230891904212e-6
            5.830003146911622e-6 -1.5685230891904212e-6 6.249736545171913
        ],
        [
            0.91889338049711 8.497430830534114 -3.2365221613872563e-6
            8.497430830534114 -3.54466336949244 -6.579625733138528e-6
            -3.2365221613872563e-6 -6.579625733138528e-6 -61.176758560765876
        ],
        [
            6.588070116109179 -0.47141097070380056 -6.329159418783092e-6
            -0.47141097070380056 -8.97451791193049 2.900790718740609e-6
            -6.329159418783092e-6 2.900790718740609e-6 61.0437074124448
        ],
        [
            -8.032517718525867 4.887120095276032 -6.45528075438051e-6
            4.887120095276032 -11.114430778391693 7.356781850376137e-6
            -6.45528075438051e-6 7.356781850376137e-6 -11.77978218702691
        ],
    ]
    @test all(isapprox.(μ, ref_μ; atol=1e-7))
end

@testitem "EffectiveMassInterpolator FourierSpaceEffectiveMass" begin
    using Wannier.Datasets
    hamiltonian, _ = read_w90_tb(dataset"Si2_valence/reference/MDRS/Si2_valence")
    interp = Wannier.EffectiveMassInterpolator(hamiltonian)

    k = [0.1, 0.2, 0.3]
    μ = interp(k, Wannier.FourierSpaceEffectiveMass())

    ref_μ = interp(k, Wannier.FiniteDifferenceEffectiveMass())
    @test all(isapprox.(μ, ref_μ; atol=2e-3))
end

@testitem "HamiltonianHessianInterpolator" begin
    using LinearAlgebra
    using Wannier.Datasets
    hamiltonian, _ = read_w90_tb(dataset"Si2_valence/reference/MDRS/Si2_valence")
    interp = Wannier.HamiltonianHessianInterpolator(hamiltonian)

    k = [0.1, 0.2, 0.3]
    d²H = interp(k)

    # only check the effmass tensor of each band, which should be real
    ref_diag = Wannier.EffectiveMassInterpolator(hamiltonian)(
        k, Wannier.FourierSpaceEffectiveMass()
    )
    # element of d²H are complex MMat3, so these also ensure that
    # the imaginary part is zero
    @test diag(d²H) ≈ ref_diag
end
