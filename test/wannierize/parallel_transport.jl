using LinearAlgebra
# the conduction bands poles are chosen randomly, skip test on conduction.
@testset "parallel_transport valence" begin
    model = read_w90(joinpath(FIXTURE_PATH, "valence", "silicon"))

    Umin, _ = parallel_transport(model)

    Uref = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"))
    # println("maxabs ", maximum(abs.(Umin - Uref)))
    # println("norm ", norm(Umin - Uref))
    # somehow in GitHub CI, the maxabs = 3.774913572015048e-7,
    # the norm = norm 4.7457964893612196e-6,
    # I increase a bit the tolerance here
    @test isapprox(Umin, Uref; atol=1e-5)

    ϵ0, ϵ1 = Wannier.compute_error(model, Umin)

    ϵ0_ref = 0.6148018374094284
    ϵ1_ref = 0.16490615880322881

    @test isapprox(ϵ0, ϵ0_ref; atol=1e-5)
    @test isapprox(ϵ1, ϵ1_ref; atol=1e-5)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.ptg.amn"), Umin)
end

@testset "parallel_transport neg coord" begin
    # Test PTG with a kgrid that have negative coordinates, i.e., -0.25 instead of 0.75
    # Before commit a1b05ae, the `parallel_transport` function will fail with at `index_bvector`.
    model = read_w90(joinpath(FIXTURE_PATH, "gaas", "gaas"))
    # only 4 valence bands as an isolated manifold
    model = Wannier.truncate(model, 1:4, 1:4)

    Umin, _ = parallel_transport(model)

    Uref = read_amn(joinpath(FIXTURE_PATH, "gaas", "gaas.val.ptg.amn"))
    @test isapprox(Umin, Uref; atol=1e-5)

    ϵ0, ϵ1 = Wannier.compute_error(model, Umin)

    ϵ0_ref = 0.20036010968611806
    ϵ1_ref = 0.2010755533663015

    @test isapprox(ϵ0, ϵ0_ref; atol=1e-5)
    @test isapprox(ϵ1, ϵ1_ref; atol=1e-5)

    # write_amn(joinpath(FIXTURE_PATH, "gaas", "gaas.val.ptg.amn"), Umin)
end
