using NLSolversBase

# A reusable fixture for a model
# no disentanglement
model = read_w90(joinpath(FIXTURE_PATH, "valence", "silicon"))
r₀ = zeros(eltype(model.lattice), 3, model.n_wann)
λ = 10.0
f, g! = Wannier.get_fg!_center_maxloc(model, r₀, λ)

@testset "constraint center maxloc spread gradient" begin
    U0 = deepcopy(model.U)

    # analytical gradient
    G = similar(U0)
    g!(G, U0)

    # finite diff gradient
    d = OnceDifferentiable(f, U0, zero(eltype(real(U0))))
    G_ref = NLSolversBase.gradient!(d, U0)

    # I am using a looser tolerance here
    @test isapprox(G, G_ref; atol=1e-6)

    # Test 2nd iteration
    # This W1 is generated from opt_rotate after 1st iteration,
    # I use it here so I don't need to store the input AMN matrices as a file.
    W1 = [
        0.497264+0.129534im -0.281405-0.540402im 0.0625545+0.237438im -0.516272-0.194678im
        0.113245-0.532338im 0.47272+0.0873106im -0.546504+0.19511im -0.36623+0.042955im
        -0.0956088+0.205024im 0.468778+0.234124im 0.563928-0.141122im -0.572565+0.0921779im
        0.561868-0.269944im 0.350754-0.00970417im 0.451343+0.247661im 0.471665-0.0282528im
    ]
    U1 = rotate_U(U0, W1)

    g!(G, U1)
    d = OnceDifferentiable(f, U1, zero(eltype(real(U1))))
    G_ref = NLSolversBase.gradient!(d, U1)
    @test isapprox(G, G_ref; atol=1e-6)
end

@testset "constraint center maxloc valence" begin
    # start from parallel transport gauge
    model.U .= read_orthonorm_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"))

    Umin = Wannier.max_localize_center(model, r₀, λ; max_iter=4)
    Ω = Wannier.omega_center(model.bvectors, model.M, Umin, r₀, λ)

    # display(Ω)
    @test Ω.Ω ≈ Ω.ΩI + Ω.Ω̃
    @test Ω.Ω̃ ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 10.79175007994904; atol=1e-7)
    @test isapprox(Ω.ΩI, 5.8127097092426245; atol=1e-7)
    @test isapprox(Ω.ΩOD, 4.912331468930677; atol=1e-7)
    @test isapprox(Ω.ΩD, 0.06670890177573785; atol=1e-7)
    @test isapprox(Ω.Ω̃, 4.979040370706415; atol=1e-7)

    @test isapprox(
        Ω.ω,
        [1.974272148108204, 2.9590704528882212, 2.957966420231007, 2.9004410587216065];
        atol=1e-7,
    )
    @test isapprox(
        Ω.r,
        [
            -0.0006371414769181565 -0.011666353819322479 0.0057784441053171055 0.00550194992572012
            0.0006811595685799503 -0.0026638301666946987 -0.0036622860182758767 0.005763601542074829
            0.0300072108960674 -0.04279486032530962 -0.03176818666470274 0.04839785532365684
        ];
        atol=1e-7,
    )

    @test Ω.Ωt ≈ Ω.Ω + Ω.Ωc
    @test isapprox(Ω.Ωc, 0.06337765901009804; atol=1e-7)
    @test isapprox(Ω.Ωt, 10.855127738959137; atol=1e-7)
    @test isapprox(
        Ω.ωc,
        [
            0.009013026333805437,
            0.01974599872857372,
            0.0105602043912133,
            0.024058429556505577,
        ];
        atol=1e-7,
    )
    @test isapprox(
        Ω.ωt,
        [1.9832851744420095, 2.9788164516167948, 2.9685266246222204, 2.924499488278112];
        atol=1e-7,
    )
    @test Ω.ωt ≈ Ω.ω + Ω.ωc
end
