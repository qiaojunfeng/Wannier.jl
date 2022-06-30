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
    d = OnceDifferentiable(f, W0, real(zero(eltype(W0))))
    G_ref = NLSolversBase.gradient!(d, W0)

    # I am using a looser tolerance here
    @test isapprox(G, G_ref; atol = 1e-6)
end


# the conduction bands poles are chosen randomly, skip test on conduction.
@testset "opt_rotate valence" begin

    # start from parallel transport gauge
    model.A .= read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"))

    Wmin = opt_rotate(model)
    Amin = rotate_amn(model.A, Wmin)

    Aref = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"))
    @test isapprox(Amin, Aref; atol = 1e-7)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.optrot.amn"), Amin)
end
