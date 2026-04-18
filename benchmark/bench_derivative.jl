module BenchDerivative

    using BenchmarkTools
    using Wannier
    using Wannier: Vec3
    using Wannier.Datasets

    SUITE = BenchmarkGroup()

    hamiltonian, _ = read_w90_tb(dataset"Si2_valence/outputs/MDRS/Si2_valence")
    # a non-degenerate kpoint
    k = [Vec3(0.1, 0.2, 0.3)]

    vel  = Wannier.VelocityInterpolator(hamiltonian)
    grad = Wannier.HamiltonianGradientInterpolator(hamiltonian)
    hess = Wannier.HamiltonianHessianInterpolator(hamiltonian)

    SUITE["velocity FD"]       = @benchmarkable $vel($k, Wannier.FiniteDifferenceVelocity())
    SUITE["velocity analytic"] = @benchmarkable $vel($k, Wannier.AnalyticVelocity())
    SUITE["hamiltonian gradient"] = @benchmarkable $grad($k)
    SUITE["hamiltonian hessian"]  = @benchmarkable $hess($k)

end  # module

BenchDerivative.SUITE
