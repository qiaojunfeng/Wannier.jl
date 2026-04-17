@testmodule SplitEnv begin
    # This code runs once and the module is cached
    using Wannier
    using Wannier.Datasets
    export model, chk, n_val, model_val, model_cond

    model = read_w90(dataset"Si2_coarse/Si2")
    chk = read_chk(dataset"Si2_coarse/outputs/Si2.chk")
    n_val = 4

    model_val = read_w90(dataset"Si2_coarse/valence/Si2_val")
    model_cond = read_w90(dataset"Si2_coarse/conduction/Si2_cond")
end

@testitem "split_eig" setup = [SplitEnv] begin
    E = model.eigenvalues
    U = WannierIO.gauge_matrices(chk)

    (Ev, _), (Ec, _) = Wannier.split_eig(E, U, n_val)

    Ev_ref = model_val.eigenvalues
    Ec_ref = model_cond.eigenvalues

    @test isapprox(Ev, Ev_ref; atol = 1.0e-7)
    @test isapprox(Ec, Ec_ref; atol = 1.0e-7)
end

@testitem "split_mmn" setup = [SplitEnv] begin
    using Wannier.Datasets

    E = model.eigenvalues
    U = WannierIO.gauge_matrices(chk)
    (Ev, _), (Ec, _) = Wannier.split_eig(E, U, n_val)

    M = model.overlaps

    # V is random, use reference V
    Vv = read_amn(dataset"Si2_coarse/valence/Si2_val.vmn").A
    Vc = read_amn(dataset"Si2_coarse/conduction/Si2_cond.vmn").A

    nk = n_kpoints(model)
    nb = n_bands(model)
    nw = n_wannier(model)
    UVv = zeros(eltype(U), nb, n_val, nk)
    UVc = zeros(eltype(U), nb, nw - n_val, nk)
    for ik in 1:nk
        view(UVv, :, :, ik) .= view(U, :, :, ik) * view(Vv, :, :, ik)
        view(UVc, :, :, ik) .= view(U, :, :, ik) * view(Vc, :, :, ik)
    end
    Mv = transform_gauge(M, kpb_k(model), UVv)
    Mc = transform_gauge(M, kpb_k(model), UVc)

    Mv_ref = model_val.overlaps
    Mc_ref = model_cond.overlaps

    @test isapprox(Mv, Mv_ref; atol = 1.0e-7)
    @test isapprox(Mc, Mc_ref; atol = 1.0e-7)

    #=
    write_eig(dataset"Si2_coarse/valence/Si2_val.eig", Ev)
    write_eig(dataset"Si2_coarse/conduction/Si2_cond.eig", Ec)
    write_amn(dataset"Si2_coarse/valence/Si2_val.vmn", Vv)
    write_amn(dataset"Si2_coarse/conduction/Si2_cond.vmn", Vc)
    write_mmn(dataset"Si2_coarse/valence/Si2_val.mmn", Mv, model.kstencil)
    write_mmn(dataset"Si2_coarse/conduction/Si2_cond.mmn", Mc, model.kstencil)
    =#
end

@testitem "identity_gauge" setup = [SplitEnv] begin
    nw = 4
    nk = 8
    U = identity_gauge(ComplexF64, nk, nw)

    Uv = model_val.gauges
    Uc = model_cond.gauges

    @test isapprox(U, Uv; atol = 1.0e-7)
    @test isapprox(U, Uc; atol = 1.0e-7)

    #=
    write_amn(dataset"Si2_coarse/valence/Si2_val.amn", U)
    write_amn(dataset"Si2_coarse/conduction/Si2_cond.amn", U)
    =#
end

@testitem "split_model" setup = [SplitEnv] begin
    model_test = deepcopy(model)
    model_test.gauges .= WannierIO.gauge_matrices(chk)

    (model_v, _), (model_c, _) = split_model(model_test, n_val)

    Ev_ref = model_val.eigenvalues
    Ec_ref = model_cond.eigenvalues

    @test isapprox(model_v.eigenvalues, Ev_ref; atol = 1.0e-7)
    @test isapprox(model_c.eigenvalues, Ec_ref; atol = 1.0e-7)
end
