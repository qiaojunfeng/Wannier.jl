import LinearAlgebra as LA
using NLSolversBase

# A reusable fixture for a model
model = read_seedname("$FIXTURE_PATH/silicon")

f, g! = Wannier.get_fg!_rotate(model)


@testset "spread gradient" begin

    W0 = LA.diagm(0 => fill(1.0 + 0 * im, model.n_wann))

    # analytical gradient
    G = similar(W0)
    g!(G, W0)

    # finite diff gradient
    d = OnceDifferentiable(f, W0, real(zero(eltype(W0))))
    G_ref = NLSolversBase.gradient!(d, W0)

    @test isapprox(G, G_ref; atol = 1e-7)
end
