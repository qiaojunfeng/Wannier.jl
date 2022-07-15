using LinearAlgebra

# the conduction bands poles are chosen randomly, skip test on conduction.
@testset "parallel_transport valence" begin
    model = read_seedname(joinpath(FIXTURE_PATH, "valence", "silicon"))

    Amin, _ = parallel_transport(model)

    Aref = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"))
    # println("maxabs ", maximum(abs.(Amin - Aref)))
    # println("norm ", norm(Amin - Aref))
    # somehow in GitHub CI, the maxabs = 3.774913572015048e-7,
    # the norm = norm 4.7457964893612196e-6,
    # I increase a bit the tolerance here
    @test isapprox(Amin, Aref; atol=1e-5)

    ϵ0, ϵ1 = Wannier.compute_error(model, Amin)

    ϵ0_ref = 0.6148018374094284
    ϵ1_ref = 0.16490615880322881

    @test isapprox(ϵ0, ϵ0_ref; atol=1e-5)
    @test isapprox(ϵ1, ϵ1_ref; atol=1e-5)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"), Amin)
end
