using NLSolversBase

# A reusable fixture for a model
# no disentanglement
model = read_w90(joinpath(FIXTURE_PATH, "valence", "silicon"))
fg! = Wannier.get_fg!_maxloc(model)

@testset "maxloc spread gradient" begin
    U0 = [model.U[ik][ib, ic] for ib=1:size(model.U[1],1), ic = 1:size(model.U[1],2), ik = 1:length(model.U)]

    # analytical gradient
    G = similar(U0)
    fg!(nothing, G, U0)

    # finite diff gradient
    d = OnceDifferentiable(x->fg!(1.0, nothing, x), U0, zero(eltype(real(U0))))
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

    fg!(nothing, G, U1)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), U1, zero(eltype(real(U1))))
    G_ref = NLSolversBase.gradient!(d, U1)
    @test isapprox(G, G_ref; atol=1e-6)
end

@testset "maxloc valence" begin
    # start from parallel transport gauge
    model.U .= read_orthonorm_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"))

    Umin = max_localize(model)
    Ω = omega(model.bvectors, model.M, Umin)

    @test isapprox(Ω.Ω, 6.374823673444644; atol=1e-7)
    @test isapprox(Ω.ΩI, 5.812709709242578; atol=1e-7)
    @test isapprox(Ω.ΩOD, 0.5621139641912166; atol=1e-7)
    @test isapprox(Ω.Ω̃, 0.5621139642020658; atol=1e-7)
end
