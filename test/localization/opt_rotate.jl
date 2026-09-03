@testmodule OptRotEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier.Datasets
    export model, fg!

    # no disentanglement
    model = read_w90(dataset"Si2_valence_coarse/Si2")
    # reset initial gauge
    model.gauges .= identity_gauge(eltype(model.gauges[1]), n_kpoints(model), n_wannier(model))
    fg! = Wannier._optimizer_callback(Wannier.Problem(Wannier.Variance(), model, Wannier.WLayout()))
end

@testitem "opt_rotate spread gradient" setup = [OptRotEnv] begin
    using LinearAlgebra
    using NLSolversBase

    W0 = diagm(0 => fill(1.0 + 0 * im, n_wannier(model)))

    # analytical gradient
    G = similar(W0)
    fg!(nothing, G, W0)

    # finite diff gradient
    d = OnceDifferentiable(W -> fg!(1.0, nothing, W), W0, zero(eltype(real(W0))))
    G_ref = NLSolversBase.gradient!(d, W0)

    # I am using a looser tolerance here
    @test isapprox(G, G_ref; atol = 1.0e-6)

    # Test 2nd iteration
    W1 = Wannier.localize(Wannier.Variance(), model, Wannier.WLayout(); max_iter = 1)

    fg!(nothing, G, W1)
    d = OnceDifferentiable(W -> fg!(1.0, nothing, W), W1, zero(eltype(real(W1))))
    G_ref = NLSolversBase.gradient!(d, W1)
    @test isapprox(G, G_ref; atol = 1.0e-6)
end

@testitem "opt_rotate valence" setup = [OptRotEnv] begin
    using Wannier.Datasets
    # reset initial gauge
    U0 = read_amn_ortho(dataset"Si2_valence_coarse/outputs/ptg.amn")
    model.gauges .= U0

    Wmin = localize(Wannier.Variance(), model, Wannier.WLayout())
    Wref = [
        0.470784 + 0.136745im    0.337519 - 0.159248im  -0.502644 - 0.281439im   -0.466234 + 0.266743im
        -0.106769 + 0.0365541im   0.755771 - 0.287907im  -0.137358 + 0.120947im      0.4286 - 0.340573im
        0.147589 - 0.213364im    0.381346 - 0.166663im   0.771383 - 0.0291013im  -0.190455 + 0.356846im
        0.605102 + 0.559885im   -0.173099 - 0.05895im    0.194089 + 0.0331043im   0.490533 + 0.0869002im
    ]
    # display(Wmin)
    @test isapprox(Wmin, Wref; atol = 1.0e-5)
end
