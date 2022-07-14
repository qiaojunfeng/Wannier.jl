import LinearAlgebra as LA
using NLSolversBase

# A reusable fixture for a model
# no disentanglement
model = read_seedname(joinpath(FIXTURE_PATH, "valence", "silicon"))
f, g! = Wannier.get_fg!_rotate(model)

@testset "opt_rotate spread gradient" begin
    W0 = LA.diagm(0 => fill(1.0 + 0 * im, model.n_wann))

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
        0.49468641+0.03617146im 0.07068715-0.51132487im 0.41874113+0.10510656im -0.47726951-0.27083108im
        0.19519914-0.37114366im 0.71429065+0.03150132im -0.41871117-0.27487794im -0.23248322+0.08955916im
        -0.63283158-0.18051594im 0.35991701+0.12405487im 0.61372044-0.01153739im -0.21265271-0.00112081im
        0.26182627-0.27660776im 0.18700395-0.20602523im 0.38307371-0.19799040im 0.76804146+0.04104790im
    ]
    @test isapprox(Wmin, Wref; atol=1e-7)

    Amin = rotate_amn(A0, Wmin)
    Aref = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"))
    @test isapprox(Amin, Aref; atol=1e-7)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"), Amin)
end
