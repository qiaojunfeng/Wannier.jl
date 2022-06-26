using YAML

@testset "read win" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/win.yaml")

    win = read_win("$FIXTURE_PATH/silicon.win")

    mat2vec(A::AbstractMatrix) = [Vector(A[:, i]) for i = 1:size(A, 2)]

    # Convert type so YAML can write it.
    win["unit_cell"] = mat2vec(win["unit_cell"])
    win["kpoints"] = mat2vec(win["kpoints"])

    # YAML.write_file(String(@__DIR__) * "/test_data/win.yaml", win)

    @test begin
        test_kpath = test_data["kpoint_path"]
        win_kpath = win["kpoint_path"]
        # println(test_kpath)
        # println(win_kpath)

        length(test_kpath) != length(win_kpath) && return false

        for i = 1:length(win_kpath)
            # a Pair: "L" => [0.5, 0.5, 0.5]
            (w1_lab, w1_vec), (w2_lab, w2_vec) = win_kpath[i]
            # YAML output is a Dict: Dict("L" => [0.5, 0.5, 0.5])
            tk1, tk2 = test_kpath[i]
            t1_lab, t1_vec = [(k, v) for (k, v) in tk1][1]
            t2_lab, t2_vec = [(k, v) for (k, v) in tk2][1]

            w1_lab != t1_lab && return false
            w1_vec ≉ t1_vec && return false
            w2_lab != t2_lab && return false
            w2_vec ≉ t2_vec && return false
        end
        true
    end

    @test test_data["num_wann"] == win["num_wann"]
    @test test_data["num_bands"] == win["num_bands"]
    @test test_data["unit_cell"] ≈ win["unit_cell"]
    @test test_data["mp_grid"] == win["mp_grid"]
    @test test_data["kpoints"] ≈ win["kpoints"]
end


@testset "read/write mmn" begin
    M, kpb_k, kpb_b = read_mmn("$FIXTURE_PATH/silicon.mmn")

    tmpfile = tempname(cleanup=true)

    write_mmn(tmpfile, M, kpb_k, kpb_b)

    M2, kpb_k2, kpb_b2 = read_mmn(tmpfile)

    @test M ≈ M2
    @test kpb_k ≈ kpb_k2
    @test kpb_b ≈ kpb_b2
end


@testset "read/write eig" begin
    E = read_eig("$FIXTURE_PATH/silicon.eig")

    tmpfile = tempname(cleanup=true)

    write_eig(tmpfile, E)

    E2 = read_eig(tmpfile)

    @test E ≈ E2
end


@testset "read_seedname" begin
    model = read_seedname("$FIXTURE_PATH/silicon")

    @test model.n_bands ≈ 12
    @test model.n_wann ≈ 8
    @test model.n_kpts ≈ 64
end
