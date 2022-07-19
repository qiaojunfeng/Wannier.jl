
@testset "split_eig" begin
    E = read_eig(joinpath(FIXTURE_PATH, "silicon/silicon.eig"))
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))
    n_val = 4

    U = get_amn(chk)

    Ev, Ec, _, _ = Wannier.split_eig(E, U, n_val)

    Ev_ref = read_eig(joinpath(FIXTURE_PATH, "valence", "silicon.eig"))
    Ec_ref = read_eig(joinpath(FIXTURE_PATH, "conduction", "silicon.eig"))

    @test isapprox(Ev, Ev_ref; atol=1e-7)
    @test isapprox(Ec, Ec_ref; atol=1e-7)
end

@testset "split_mmn" begin
    E = read_eig(joinpath(FIXTURE_PATH, "silicon/silicon.eig"))
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))
    n_val = 4

    U = get_amn(chk)
    Ev, Ec, _, _ = Wannier.split_eig(E, U, n_val)

    M, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "silicon/silicon.mmn"))

    # V is random, use reference V
    Vv = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.vmn"))
    Vc = read_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.vmn"))

    n_kpts = size(M, 4)
    n_bands = size(M, 1)
    n_wann = size(U, 2)
    UVv = similar(U, n_bands, n_val, n_kpts)
    UVc = similar(U, n_bands, n_wann - n_val, n_kpts)
    for ik in 1:n_kpts
        UVv[:, :, ik] = U[:, :, ik] * Vv[:, :, ik]
        UVc[:, :, ik] = U[:, :, ik] * Vc[:, :, ik]
    end
    Mv = rotate_M(M, kpb_k, UVv)
    Mc = rotate_M(M, kpb_k, UVc)

    Mv_ref, _, _ = read_mmn(joinpath(FIXTURE_PATH, "valence", "silicon.mmn"))
    Mc_ref, _, _ = read_mmn(joinpath(FIXTURE_PATH, "conduction", "silicon.mmn"))

    @test isapprox(Mv, Mv_ref; atol=1e-7)
    @test isapprox(Mc, Mc_ref; atol=1e-7)

    # write_eig(joinpath(FIXTURE_PATH, "valence", "silicon.eig"), Ev)
    # write_eig(joinpath(FIXTURE_PATH, "conduction", "silicon.eig"), Ec)
    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.vmn"), Vv)
    # write_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.vmn"), Vc)
    # write_mmn(joinpath(FIXTURE_PATH, "valence", "silicon.mmn"), Mv, kpb_k, kpb_b)
    # write_mmn(joinpath(FIXTURE_PATH, "conduction", "silicon.mmn"), Mc, kpb_k, kpb_b)
end

@testset "eyes_A" begin
    n_wann = 4
    n_kpts = 64

    A = eyes_A(ComplexF64, n_wann, n_kpts)

    Av = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.amn"))
    Ac = read_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.amn"))

    @test isapprox(A, Av; atol=1e-7)
    @test isapprox(A, Ac; atol=1e-7)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.amn"), A)
    # write_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.amn"), A)
end

@testset "split_model" begin
    model = read_seedname(joinpath(FIXTURE_PATH, "silicon/silicon"))
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))
    model.A .= get_amn(chk)

    n_val = 4

    model_v, model_c = split_model(model, n_val)

    Ev_ref = read_eig(joinpath(FIXTURE_PATH, "valence", "silicon.eig"))
    Ec_ref = read_eig(joinpath(FIXTURE_PATH, "conduction", "silicon.eig"))

    @test isapprox(model_v.E, Ev_ref; atol=1e-7)
    @test isapprox(model_c.E, Ec_ref; atol=1e-7)
end
