using LinearAlgebra
using NLSolversBase

# A reusable fixture for a model
# no disentanglement
model = read_w90(joinpath(FIXTURE_PATH, "valence", "silicon"))
f, g! = Wannier.get_fg!_rotate(model)

@testset "opt_rotate spread gradient" begin
    W0 = diagm(0 => fill(1.0 + 0 * im, model.n_wann))

    # analytical gradient
    G = similar(W0)
    g!(G, W0)

    # finite diff gradient
    d = OnceDifferentiable(f, W0, zero(eltype(real(W0))))
    G_ref = NLSolversBase.gradient!(d, W0)

    # I am using a looser tolerance here
    @test isapprox(G, G_ref; atol=1e-6)

    # Test 2nd iteration
    W1 = [
        0.497264+0.129534im -0.281405-0.540402im 0.0625545+0.237438im -0.516272-0.194678im
        0.113245-0.532338im 0.47272+0.0873106im -0.546504+0.19511im -0.36623+0.042955im
        -0.0956088+0.205024im 0.468778+0.234124im 0.563928-0.141122im -0.572565+0.0921779im
        0.561868-0.269944im 0.350754-0.00970417im 0.451343+0.247661im 0.471665-0.0282528im
    ]

    g!(G, W1)
    d = OnceDifferentiable(f, W1, zero(eltype(real(W1))))
    G_ref = NLSolversBase.gradient!(d, W1)
    @test isapprox(G, G_ref; atol=1e-6)
end

@testset "opt_rotate valence" begin
    # start from parallel transport gauge
    A0 = read_orthonorm_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"))
    model.A .= A0

    Wmin = opt_rotate(model)
    Wref = [
        0.494686+0.036171im 0.070687-0.511325im 0.418741+0.105107im -0.477270-0.270831im
        0.195199-0.371144im 0.714291+0.031501im -0.418711-0.274878im -0.232483+0.089559im
        -0.632832-0.180516im 0.359917+0.124055im 0.613720-0.011537im -0.212653-0.001121im
        0.261826-0.276608im 0.187004-0.206025im 0.383074-0.197990im 0.768041+0.041048im
    ]
    # display(Wmin)
    @test isapprox(Wmin, Wref; atol=1e-5)

    # Amin = rotate_A(A0, Wmin)
    # Aref = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"))
    # @test isapprox(Amin, Aref; atol=1e-6)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"), Amin)
end
