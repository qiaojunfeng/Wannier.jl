import LinearAlgebra as LA


@testset "split_eig" begin
    E = read_eig(joinpath(FIXTURE_PATH, "silicon.eig"))
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon.chk.fmt"))
    n_val = 4

    Ev, Ec, U, V = Wannier.split_eig(E, chk, n_val)

    Ev_ref = read_eig(joinpath(FIXTURE_PATH, "valence", "silicon.eig"))
    Ec_ref = read_eig(joinpath(FIXTURE_PATH, "conduction", "silicon.eig"))

    @test isapprox(Ev, Ev_ref; atol = 1e-7)
    @test isapprox(Ec, Ec_ref; atol = 1e-7)
end


@testset "split_mmn" begin
    E = read_eig(joinpath(FIXTURE_PATH, "silicon.eig"))
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon.chk.fmt"))
    n_val = 4

    Ev, Ec, U, V = Wannier.split_eig(E, chk, n_val)

    M, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "silicon.mmn"))

    # V is random, use reference V
    Vv = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.vmn"))
    Vc = read_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.vmn"))
    V_ref = cat(Vv, Vc, dims = 2)

    Mv, Mc = Wannier.split_mmn(M, kpb_k, U, V_ref, n_val)

    Mv_ref, _, _ = read_mmn(joinpath(FIXTURE_PATH, "valence", "silicon.mmn"))
    Mc_ref, _, _ = read_mmn(joinpath(FIXTURE_PATH, "conduction", "silicon.mmn"))

    @test isapprox(Mv, Mv_ref; atol = 1e-7)
    @test isapprox(Mc, Mc_ref; atol = 1e-7)

    # write_eig(joinpath(FIXTURE_PATH, "valence", "silicon.eig"), Ev)
    # write_eig(joinpath(FIXTURE_PATH, "conduction", "silicon.eig"), Ec)
    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.vmn"), Vv)
    # write_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.vmn"), Vc)
    # write_mmn(joinpath(FIXTURE_PATH, "valence", "silicon.mmn"), Mv, kpb_k, kpb_b)
    # write_mmn(joinpath(FIXTURE_PATH, "conduction", "silicon.mmn"), Mc, kpb_k, kpb_b)
end


@testset "ones_amn" begin
    n_wann = 4
    n_kpts = 64

    A = ones_amn(ComplexF64, n_wann, n_kpts)

    Av = read_amn(joinpath(FIXTURE_PATH, "valence", "silicon.amn"))
    Ac = read_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.amn"))

    @test isapprox(A, Av; atol = 1e-7)
    @test isapprox(A, Ac; atol = 1e-7)

    # write_amn(joinpath(FIXTURE_PATH, "valence", "silicon.amn"), A)
    # write_amn(joinpath(FIXTURE_PATH, "conduction", "silicon.amn"), A)
end


@testset "split_valence_conduction" begin
    tmpdir_val = mktempdir(cleanup = true)
    tmpdir_cond = mktempdir(cleanup = true)

    n_val = 4
    split_valence_conduction(
        joinpath(FIXTURE_PATH, "silicon"),
        n_val,
        unk = false,
        vmn = true,
        outdir_val = tmpdir_val,
        outdir_cond = tmpdir_cond,
    )

    Ev = read_eig(joinpath(tmpdir_val, "silicon.eig"))
    Ec = read_eig(joinpath(tmpdir_cond, "silicon.eig"))

    Ev_ref = read_eig(joinpath(FIXTURE_PATH, "valence", "silicon.eig"))
    Ec_ref = read_eig(joinpath(FIXTURE_PATH, "conduction", "silicon.eig"))

    @test isapprox(Ev, Ev_ref; atol = 1e-7)
    @test isapprox(Ec, Ec_ref; atol = 1e-7)
end
