@testitem "parity: Variance disentangle matches legacy path" begin
    using Wannier
    using Wannier.Datasets

    model = read_w90(dataset"Si2_coarse/Si2")

    # Legacy path
    U_legacy = Wannier.disentangle(model; max_iter = 10, f_tol = 1.0e-8, g_tol = 1.0e-6)

    # New path: Problem + solve!
    prob = Wannier.Problem(Wannier.Variance(), model)
    solver = Wannier.OptimLBFGS(; max_iter = 10, f_tol = 1.0e-8, g_tol = 1.0e-6)
    U_new = Wannier.solve!(prob, solver)

    # Delegated fg! ⇒ byte-identical trajectory
    @test isapprox(U_legacy, U_new; atol = 1.0e-12)
end

@testitem "parity: Variance max_localize matches legacy path" begin
    using Wannier
    using Wannier.Datasets

    model = read_w90_with_chk(dataset"Si2_valence_coarse/Si2", dataset"Si2_valence_coarse/outputs/Si2.chk")

    U_legacy = Wannier.max_localize(model; max_iter = 10, f_tol = 1.0e-8, g_tol = 1.0e-6)

    prob = Wannier.Problem(Wannier.Variance(), model)
    solver = Wannier.OptimLBFGS(; max_iter = 10, f_tol = 1.0e-8, g_tol = 1.0e-6)
    U_new = Wannier.solve!(prob, solver)

    @test isapprox(U_legacy, U_new; atol = 1.0e-12)
end
