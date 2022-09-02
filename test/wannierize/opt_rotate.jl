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
        0.49473+0.0355699im 0.0673789-0.511771im 0.418917+0.104404im -0.478045-0.26946im
        0.194747-0.371381im 0.71448+0.0268801im -0.419172-0.274175im -0.232225+0.090227im
        -0.633051-0.179746im 0.360712+0.121724im 0.6137-0.012567im -0.212654-0.000510533im
        0.26149-0.276926im 0.185667-0.207231im 0.382741-0.198633im 0.768156+0.0388423im
    ]
    # display(Wmin)
    @test isapprox(Wmin, Wref; atol=1e-5)

    # Amin = rotate_A(A0, Wmin)
    # Aref = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"))
    # @test isapprox(Amin, Aref; atol=1e-6)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"), Amin)
end
