@testmodule DisentangleEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier.Datasets
    export model, fg!

    model = read_w90(dataset"Si2_coarse/Si2")
    fg! = Wannier._optimizer_callback(Wannier.Problem(Wannier.Variance(), model, Wannier.XYLayout()))
end

@testitem "compact XY layout" setup = [DisentangleEnv] begin
    layout = Wannier.XYLayout()
    @test Wannier.Problem(Wannier.Variance(), model, layout).evaluation === nothing
    x = Wannier.initial_parameters(layout, model)
    U0 = Wannier.finalize_result(layout, x, model)

    # Initialization is stable at the canonical-gauge interface even though
    # the internal X/Y factorization is nonunique.
    model1 = deepcopy(model)
    model1.gauges .= U0
    x1 = Wannier.initial_parameters(layout, model1)
    U1 = Wannier.finalize_result(layout, x1, model1)
    @test isapprox(U0, U1; atol = 1.0e-6)

    expected = sum(1:n_kpoints(model)) do ik
        nf = count(view(model.frozen_bands, :, ik))
        n_wannier(model)^2 + (n_bands(model) - nf) * (n_wannier(model) - nf)
    end
    @test length(x) == expected
end

@testitem "compact XY with arbitrary frozen rows" begin
    using LinearAlgebra
    using Random

    rng = MersenneTwister(1234)
    nbands, nwann, nkpts = 6, 3, 4
    frozen = falses(nbands, nkpts)
    frozen[2, 2] = true
    frozen[[1, 3, 5], 3] .= true
    frozen[1:3, 4] .= true
    xy = Wannier._xy_structure(frozen, nwann)

    x = zeros(ComplexF64, xy.nparameters)
    for block in xy.blocks
        X = reshape(view(x, block.x_range), nwann, nwann)
        X .= Matrix(qr(randn(rng, ComplexF64, nwann, nwann)).Q)
        nfrozen = length(block.frozen)
        nactive = nwann - nfrozen
        if nactive > 0
            Y = reshape(view(x, block.y_range), length(block.nonfrozen), nactive)
            Y .= Matrix(
                qr(randn(rng, ComplexF64, length(block.nonfrozen), nactive)).Q
            )[:, 1:nactive]
        end
    end

    X = zeros(ComplexF64, nwann, nwann, nkpts)
    Y = zeros(ComplexF64, nbands, nwann, nkpts)
    U = similar(Y)
    Wannier._initialize_compact_y!(Y, xy)
    Wannier.assemble_gauge!(U, X, Y, x, xy)
    @test_throws DimensionMismatch Wannier.assemble_gauge!(U, X, Y, x[1:(end - 1)], xy)

    for (ik, block) in enumerate(xy.blocks)
        nfrozen = length(block.frozen)
        nactive = nwann - nfrozen
        Xk = reshape(view(x, block.x_range), nwann, nwann)
        @test U[block.frozen, :, ik] ≈ Xk[1:nfrozen, :]
        if nactive > 0
            Yk = reshape(view(x, block.y_range), length(block.nonfrozen), nactive)
            @test U[block.nonfrozen, :, ik] ≈ Yk * Xk[(nfrozen + 1):nwann, :]
        end
    end

    GU = randn(rng, ComplexF64, size(U))
    g = similar(x)
    Wannier.pullback_gradient!(g, GU, X, Y, xy)
    dx = randn(rng, ComplexF64, length(x))
    dx ./= norm(dx)
    function linear_form(z)
        Xz = similar(X)
        Yz = similar(Y)
        Uz = similar(U)
        Wannier._initialize_compact_y!(Yz, xy)
        Wannier.assemble_gauge!(Uz, Xz, Yz, z, xy)
        return 2 * real(sum(conj.(GU) .* Uz))
    end
    ε = 1.0e-6
    fd = (linear_form(x .+ ε .* dx) - linear_form(x .- ε .* dx)) / (2ε)
    analytical = 2 * real(sum(conj.(g) .* dx))
    @test isapprox(fd, analytical; rtol = 1.0e-8)

    manifold = Wannier.XYManifold(xy)
    xtrial = x + 0.1 * randn(rng, ComplexF64, length(x))
    Wannier.Optim.retract!(manifold, xtrial)
    for block in manifold.blocks
        Xk = reshape(view(xtrial, block.x_range), block.nwann, block.nwann)
        @test Xk' * Xk ≈ I
        if !isempty(block.y_range)
            Yk = reshape(view(xtrial, block.y_range), block.ynrows, block.yncols)
            @test Yk' * Yk ≈ I
        end
    end
end

@testitem "disentangle spread gradient" setup = [DisentangleEnv] begin
    using NLSolversBase

    # analytical gradient
    layout = Wannier.XYLayout()
    XY = Wannier.initial_parameters(layout, model)
    G = similar(XY)
    fg!(nothing, G, XY)

    # finite diff gradient
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)

    @test isapprox(G, G_ref; atol = 1.0e-6)

    # Test 2nd iteration
    U1 = Wannier.localize(model; max_iter = 1)
    model1 = deepcopy(model)
    model1.gauges .= U1
    XY = Wannier.initial_parameters(layout, model1)

    fg!(nothing, G, XY)
    d = OnceDifferentiable(x -> fg!(1.0, nothing, x), XY, zero(eltype(real(XY))))
    G_ref = NLSolversBase.gradient!(d, XY)
    @test isapprox(G, G_ref; atol = 1.0e-6)
end

@testitem "disentangle" setup = [DisentangleEnv] begin
    Umin = Wannier.localize(model; max_iter = 4)
    Ω = Wannier.spread(model.kstencil, model.overlaps, Umin)

    # display(Ω)
    @test Ω.Ω ≈ Ω.ΩI + Ω.Ωtilde
    @test Ω.Ωtilde ≈ Ω.ΩOD + Ω.ΩD
    @test isapprox(Ω.Ω, 12.362335109447647; atol = 1.0e-7)
    @test isapprox(Ω.ΩI, 7.212573765139664; atol = 1.0e-7)
    @test isapprox(Ω.ΩOD, 4.929849448152594; atol = 1.0e-7)
    @test isapprox(Ω.ΩD, 0.2199118961553884; atol = 1.0e-7)

    @test isapprox(
        Ω.ω,
        [
            1.2145730431631943,
            1.655531487632798,
            1.655531511083343,
            1.6555315112860653,
            1.2145730436877367,
            1.6555314861440165,
            1.6555315134652924,
            1.6555315129852053,
        ];
        atol = 1.0e-7,
    )
    @test isapprox(
        Ω.r,
        [
            [9.386119783477958e-9, -4.499475172965251e-8, 2.4981053542113326e-8],
            [-4.567961841769885e-8, 1.3926293592272931e-8, -4.076140836146991e-8],
            [-6.90776277038549e-8, -3.595155572896913e-8, -5.210872698786784e-8],
            [6.506912168541962e-8, 1.564160769284076e-7, 3.341713422557596e-8],
            [1.3576325402117597, 1.3576324575774867, 1.3576324919985168],
            [1.3576325378090364, 1.357632474151124, 1.357632542136112],
            [1.357632344302077, 1.3576324993020412, 1.3576325171624346],
            [1.3576324580526173, 1.3576326584763287, 1.3576324408622935],
        ];
        atol = 1.0e-7,
    )
end
