using Wannier: Vec3
model = read_w90_tb(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))
Rvecs = model.R.Rvectors
H = model.H
# choose a random k-point such that there is no degeneracy
k = [Vec3(0.1, 0.2, 0.3)]

@testset "velocity_fd" begin
    V = Wannier.velocity_fd(Rvecs, H, k)  # n_wann * n_kpts * 3

    # n_wann * 3, kpoint index is ignored since it is 1
    ref_V = [
        -1.1626962147528452 3.185414884939508 0.5133590506791741
        4.821069116725196 -5.903791392503299 0.18504317131351478
        3.46519032366599 -5.257640432705379 -0.2905786249951525
        -1.9838751726628523 -4.521763546720958 -0.4078235615945225
    ]

    @test all(isapprox.(V[:, 1, :], ref_V; atol=1e-7))
end

@testset "velocity" begin
    E, V = Wannier.velocity(H, k)

    ref_E = [
        -4.9062394299442555
        3.0236804371409733
        4.52689293393807
        5.5325326771182715
    ]
    @test all(isapprox.(E[1], ref_E; atol=1e-7))

    # n_wann * 3
    ref_V = [
        -1.162691011109657 3.185415086017154 0.513361677372566
        4.821074793026908 -5.903788558835394 0.18504928037519264
        3.4651970000631547 -5.257651521323197 -0.2905853774332443
        -1.9838856106328957 -4.521775180868606 -0.4078255449019925
    ]

    @test all(isapprox.(V[:, 1, :], ref_V; atol=1e-7))
end

@testset "get_dH_da" begin
    V = Wannier.get_dH_da(H, k)  # n_wann * n_wann * 1 * 3

    # n_wann * 3
    ref_Vdiag = [
        -1.162691011109657 3.185415086017154 0.513361677372566
        4.821074793026908 -5.903788558835394 0.18504928037519264
        3.4651970000631547 -5.257651521323197 -0.2905853774332443
        -1.9838856106328957 -4.521775180868606 -0.4078255449019925
    ]

    Vx = V[:, :, 1, 1]
    Vy = V[:, :, 1, 2]
    Vz = V[:, :, 1, 3]
    @test all(isapprox.(diag(Vx), ref_Vdiag[:, 1]; atol=1e-7))
    @test all(isapprox.(diag(Vy), ref_Vdiag[:, 2]; atol=1e-7))
    @test all(isapprox.(diag(Vz), ref_Vdiag[:, 3]; atol=1e-7))

    Vx_offdiag = Vx - diagm(diag(Vx))
    Vy_offdiag = Vy - diagm(diag(Vy))
    Vz_offdiag = Vz - diagm(diag(Vz))
    @test all(isapprox.(Vx_offdiag, 0; atol=1e-10))
    @test all(isapprox.(Vy_offdiag, 0; atol=1e-10))
    @test all(isapprox.(Vz_offdiag, 0; atol=1e-10))
end

@testset "effmass_fd" begin
    μ = Wannier.effmass_fd(model.R, H, k)  # n_wann * n_kpts * 3 * 3

    ref_μ = zeros(size(model.H[1].block, 1), 3, 3)
    ref_μ[1, :, :] = [
        4.554781130928 2.556346799132 -1.537850125466
        2.556346799132 5.257344720455 0.699701995011
        -1.537850125466 0.699701995011 6.283941746510
    ]
    ref_μ[2, :, :] = [
        -0.344714405909 -9.746684513701 0.834462826838
        -9.746684513701 4.755988562533 3.076987332484
        0.834462826838 3.076987332484 -54.592262209052
    ]
    ref_μ[3, :, :] = [
        -7.631248695184 -0.709917980402 1.382780718373
        -0.709917980402 1.695729110907 -3.550658202123
        1.382780718373 -3.550658202123 42.420041201474
    ]
    ref_μ[4, :, :] = [
        -14.936209280059 -7.980357962900 -0.679393347802
        -7.980357962900 -10.229698156117 -0.226031142248
        -0.679393347802 -0.226031142248 -7.805563588370
    ]

    @test all(isapprox.(μ[:, 1, :, :], ref_μ; atol=1e-7))
end

@testset "get_d2H_dadb" begin
    V = Wannier.get_d2H_dadb(H, k)  # n_wann * n_wann * 1 * 3 * 3

    # n_wann * 3 * 3, I only check the effmass tensor of each band, which should be real.
    ref_V = zeros(Wannier.n_wann(model.H), 3, 3)
    ref_V[1, :, :] = [
        4.554769651183 2.556355324663 -1.537854851936
        2.556355324663 5.257336684714 0.699697328906
        -1.537854851936 0.699697328906 6.283937437599
    ]
    ref_V[2, :, :] = [
        -0.344655021796 -9.746693370131 0.834467149282
        -9.746693370131 4.756025559528 3.077025103346
        0.834467149282 3.077025103346 -54.59320163663
    ]
    ref_V[3, :, :] = [
        -7.631243789375 -0.709921392767 1.382778765576
        -0.709921392767 1.695737800221 -3.550695766566
        1.382778765576 -3.550695766566 42.421002608166
    ]
    ref_V[4, :, :] = [
        -14.936269464782 -7.980376014297 -0.679390997763
        -7.980376014297 -10.229715226381 -0.226026667223
        -0.679390997763 -0.226026667223 -7.805582561336
    ]

    V1 = V[1, 1, 1, :, :]
    V2 = V[2, 2, 1, :, :]
    V3 = V[3, 3, 1, :, :]
    V4 = V[4, 4, 1, :, :]
    # V1, ..., V4 are complex matrices, so these also ensure that
    # the imaginary part is zero.
    @test all(isapprox.(V1, ref_V[1, :, :]; atol=1e-7))
    @test all(isapprox.(V2, ref_V[2, :, :]; atol=1e-7))
    @test all(isapprox.(V3, ref_V[3, :, :]; atol=1e-7))
    @test all(isapprox.(V4, ref_V[4, :, :]; atol=1e-7))
end
